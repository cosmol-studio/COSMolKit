//! Shared topology value container.
//!
//! This type stores concrete graph rows and local adjacency/annotation data.
//! Topology remapping policies and live-molecule installation remain owned by
//! the runtime crate.

use crate::{
    AdjacencyList, Atom, AtomId, AtomMapping, AtomSpec, Bond, BondId, BondMapping, BondSpec,
    BondStereo, BondValueError, MappingValidationError, StereoGroup, SubstanceGroup,
    TemplateAttachmentOrderError, TopologyMapping,
};

#[derive(Debug, Clone, PartialEq, Eq, thiserror::Error)]
pub enum TopologyValidationError {
    #[error("atom at position {position} has id {id}, expected {position}")]
    AtomIdMismatch { position: usize, id: AtomId },
    #[error("atom {atom} has invalid template attachment order: {source}")]
    TemplateAttachmentOrder {
        atom: AtomId,
        source: TemplateAttachmentOrderError,
    },
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

#[derive(Debug, Clone, PartialEq, Eq, thiserror::Error)]
pub enum TopologyEditError {
    #[error("invalid source topology: {0}")]
    InvalidSource(TopologyValidationError),
    #[error("atom {atom} is out of range for {atom_count} atoms")]
    AtomOutOfRange { atom: AtomId, atom_count: usize },
    #[error("bond {bond} is out of range for {bond_count} bonds")]
    BondOutOfRange { bond: BondId, bond_count: usize },
    #[error("bond {begin}-{end} already exists")]
    DuplicateBond { begin: AtomId, end: AtomId },
    #[error("invalid bond: {0}")]
    InvalidBond(BondValueError),
    #[error("atom permutation has length {actual}, expected {expected}")]
    PermutationLength { actual: usize, expected: usize },
    #[error(
        "atom permutation position {position} references atom {atom}, out of range for {atom_count} atoms"
    )]
    PermutationAtomOutOfRange {
        position: usize,
        atom: AtomId,
        atom_count: usize,
    },
    #[error("atom permutation position {position} repeats atom {atom}")]
    PermutationDuplicateAtom { position: usize, atom: AtomId },
    #[error("invalid topology mapping: {0}")]
    InvalidMapping(MappingValidationError),
    #[error("atom {carrier} template attachment remap failed: {source}")]
    TemplateAttachmentRemap {
        carrier: AtomId,
        source: TemplateAttachmentOrderError,
    },
    #[error("invalid edited topology: {0}")]
    InvalidResult(TopologyValidationError),
}

/// Owned edit state for a detached topology value.
///
/// The source copy is retained only to establish the mapping and abort/failure
/// isolation. This type has no authority over a live `Molecule`.
#[derive(Debug, Clone)]
pub struct TopologyBatchEdit {
    source: TopologyBlock,
    working: TopologyBlock,
    remove_atoms: Vec<bool>,
    remove_bonds: Vec<bool>,
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
    pub fn try_from_parts(
        atoms: Vec<Atom>,
        bonds: Vec<Bond>,
        substance_groups: Vec<SubstanceGroup>,
        stereo_groups: Vec<StereoGroup>,
    ) -> Result<Self, TopologyValidationError> {
        let adjacency =
            AdjacencyList::try_from_topology(atoms.len(), &bonds).map_err(|error| match error {
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
                crate::AdjacencyError::DuplicateBondId { .. }
                | crate::AdjacencyError::DuplicateEdge { .. } => {
                    TopologyValidationError::AdjacencyMismatch
                }
            })?;
        let topology = Self {
            atoms,
            bonds,
            adjacency,
            substance_groups,
            stereo_groups,
        };
        topology.validate()?;
        Ok(topology)
    }

    pub fn begin_batch_edit(&self) -> Result<TopologyBatchEdit, TopologyEditError> {
        // RDKit✔️❌: dp_delAtoms.reset(new boost::dynamic_bitset<>(getNumAtoms()));
        // RDKit✔️❌: dp_delBonds.reset(new boost::dynamic_bitset<>(getNumBonds()));
        //
        // The detached editor deliberately clones the source so failure and
        // abort cannot publish partial state. This is O(n) instead of RDKit's
        // in-place O(n) bitset setup, hence the performance marker.
        self.validate().map_err(TopologyEditError::InvalidSource)?;
        Ok(TopologyBatchEdit {
            source: self.clone(),
            working: self.clone(),
            remove_atoms: vec![false; self.atoms.len()],
            remove_bonds: vec![false; self.bonds.len()],
        })
    }

    pub fn reordered_atoms(
        &self,
        old_atom_order: &[AtomId],
    ) -> Result<(Self, TopologyMapping), TopologyEditError> {
        self.validate().map_err(TopologyEditError::InvalidSource)?;
        if old_atom_order.len() != self.atoms.len() {
            return Err(TopologyEditError::PermutationLength {
                actual: old_atom_order.len(),
                expected: self.atoms.len(),
            });
        }
        let mut seen = vec![false; self.atoms.len()];
        for (position, atom) in old_atom_order.iter().copied().enumerate() {
            if atom.index() >= self.atoms.len() {
                return Err(TopologyEditError::PermutationAtomOutOfRange {
                    position,
                    atom,
                    atom_count: self.atoms.len(),
                });
            }
            if seen[atom.index()] {
                return Err(TopologyEditError::PermutationDuplicateAtom { position, atom });
            }
            seen[atom.index()] = true;
        }

        let mut atom_old_to_new = vec![None; self.atoms.len()];
        for (new_index, old_id) in old_atom_order.iter().copied().enumerate() {
            let new_id = AtomId::new(new_index);
            atom_old_to_new[old_id.index()] = Some(new_id);
        }
        let mut atoms = Vec::with_capacity(self.atoms.len());
        for (new_index, old_id) in old_atom_order.iter().copied().enumerate() {
            let mut atom = self.atoms[old_id.index()]
                .clone()
                .with_id(AtomId::new(new_index));
            atom.remap_template_attachment_order(&atom_old_to_new)
                .map_err(|source| TopologyEditError::TemplateAttachmentRemap {
                    carrier: old_id,
                    source,
                })?;
            atoms.push(atom);
        }
        let atom_new_to_old = old_atom_order.iter().copied().map(Some).collect();
        let bond_old_to_new = (0..self.bonds.len())
            .map(|index| Some(BondId::new(index)))
            .collect::<Vec<_>>();
        let bond_new_to_old = bond_old_to_new.clone();
        let mut bonds = Vec::with_capacity(self.bonds.len());
        for bond in &self.bonds {
            let begin = atom_old_to_new
                .get(bond.begin().index())
                .and_then(|mapped| *mapped)
                .ok_or_else(|| {
                    TopologyEditError::InvalidSource(
                        TopologyValidationError::BondEndpointOutOfRange {
                            bond: bond.id(),
                            endpoint: "begin",
                            atom: bond.begin(),
                            atom_count: self.atoms.len(),
                        },
                    )
                })?;
            let end = atom_old_to_new
                .get(bond.end().index())
                .and_then(|mapped| *mapped)
                .ok_or_else(|| {
                    TopologyEditError::InvalidSource(
                        TopologyValidationError::BondEndpointOutOfRange {
                            bond: bond.id(),
                            endpoint: "end",
                            atom: bond.end(),
                            atom_count: self.atoms.len(),
                        },
                    )
                })?;
            let stereo_atoms = match bond.stereo_atoms() {
                Some([left, right]) => Some([
                    atom_old_to_new
                        .get(left.index())
                        .and_then(|mapped| *mapped)
                        .ok_or_else(|| {
                            TopologyEditError::InvalidSource(
                                TopologyValidationError::StereoAtomOutOfRange {
                                    bond: bond.id(),
                                    begin: left,
                                    end: right,
                                    atom_count: self.atoms.len(),
                                },
                            )
                        })?,
                    atom_old_to_new
                        .get(right.index())
                        .and_then(|mapped| *mapped)
                        .ok_or_else(|| {
                            TopologyEditError::InvalidSource(
                                TopologyValidationError::StereoAtomOutOfRange {
                                    bond: bond.id(),
                                    begin: left,
                                    end: right,
                                    atom_count: self.atoms.len(),
                                },
                            )
                        })?,
                ]),
                None => None,
            };
            bonds.push(bond.clone().remapped(bond.id(), begin, end, stereo_atoms));
        }
        let sgroup_map = (0..self.substance_groups.len())
            .map(crate::SubstanceGroupId::new)
            .map(Some)
            .collect::<Vec<_>>();
        let substance_groups = self
            .substance_groups
            .iter()
            .enumerate()
            .filter_map(|(index, group)| {
                group.remapped(
                    crate::SubstanceGroupId::new(index),
                    &atom_old_to_new,
                    &bond_old_to_new,
                    &sgroup_map,
                )
            })
            .collect();
        let stereo_groups = self
            .stereo_groups
            .iter()
            .filter_map(|group| group.remapped(&atom_old_to_new, &bond_old_to_new))
            .collect();
        let mapping = TopologyMapping {
            atoms: AtomMapping {
                old_to_new: atom_old_to_new,
                new_to_old: atom_new_to_old,
            },
            bonds: BondMapping {
                old_to_new: bond_old_to_new,
                new_to_old: bond_new_to_old,
            },
        };
        mapping
            .validate_for_counts(self.atoms.len(), atoms.len(), self.bonds.len(), bonds.len())
            .map_err(TopologyEditError::InvalidMapping)?;
        let topology = Self::try_from_parts(atoms, bonds, substance_groups, stereo_groups)
            .map_err(TopologyEditError::InvalidResult)?;
        Ok((topology, mapping))
    }

    pub fn validate(&self) -> Result<(), TopologyValidationError> {
        for (position, atom) in self.atoms.iter().enumerate() {
            if atom.id() != AtomId::new(position) {
                return Err(TopologyValidationError::AtomIdMismatch {
                    position,
                    id: atom.id(),
                });
            }
            if let Some(order) = atom.template_attachment_order() {
                order
                    .validate_for_atom_count(self.atoms.len())
                    .map_err(|source| TopologyValidationError::TemplateAttachmentOrder {
                        atom: atom.id(),
                        source,
                    })?;
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
                    crate::AdjacencyError::DuplicateBondId { .. }
                    | crate::AdjacencyError::DuplicateEdge { .. } => {
                        TopologyValidationError::AdjacencyMismatch
                    }
                }
            })?;
        if self.adjacency != expected {
            return Err(TopologyValidationError::AdjacencyMismatch);
        }
        Ok(())
    }
}

impl TopologyBatchEdit {
    pub fn add_atom(&mut self, spec: AtomSpec) -> AtomId {
        // RDKit✔️✔️: if (dp_delAtoms->size() < getNumAtoms()) {
        // RDKit✔️✔️:   dp_delAtoms->resize(getNumAtoms());
        // RDKit✔️✔️: }
        let id = AtomId::new(self.working.atoms.len());
        self.working.atoms.push(Atom::from_spec(id, spec));
        self.remove_atoms.push(false);
        id
    }

    pub fn add_bond(&mut self, spec: BondSpec) -> Result<BondId, TopologyEditError> {
        spec.validate().map_err(TopologyEditError::InvalidBond)?;
        for atom in [spec.begin(), spec.end()] {
            if atom.index() >= self.working.atoms.len() {
                return Err(TopologyEditError::AtomOutOfRange {
                    atom,
                    atom_count: self.working.atoms.len(),
                });
            }
        }
        if let Some([begin, end]) = spec.stereo_atoms()
            && (begin.index() >= self.working.atoms.len()
                || end.index() >= self.working.atoms.len())
        {
            return Err(TopologyEditError::InvalidResult(
                TopologyValidationError::StereoAtomOutOfRange {
                    bond: BondId::new(self.working.bonds.len()),
                    begin,
                    end,
                    atom_count: self.working.atoms.len(),
                },
            ));
        }
        if spec.begin() == spec.end() {
            return Err(TopologyEditError::InvalidResult(
                TopologyValidationError::SelfLoopBond {
                    bond: BondId::new(self.working.bonds.len()),
                    atom: spec.begin(),
                },
            ));
        }
        if self.working.bonds.iter().any(|bond| {
            (bond.begin() == spec.begin() && bond.end() == spec.end())
                || (bond.begin() == spec.end() && bond.end() == spec.begin())
        }) {
            return Err(TopologyEditError::DuplicateBond {
                begin: spec.begin(),
                end: spec.end(),
            });
        }
        let id = BondId::new(self.working.bonds.len());
        self.working.bonds.push(Bond::from_spec(id, spec));
        self.remove_bonds.push(false);
        Ok(id)
    }

    pub fn remove_atom(&mut self, atom: AtomId) -> Result<(), TopologyEditError> {
        // RDKit✔️✔️: void RWMol::removeAtom(unsigned int idx) {
        // RDKit✔️✔️:   removeAtom(getAtomWithIdx(idx));
        // RDKit✔️✔️: }
        let atom_count = self.working.atoms.len();
        let Some(slot) = self.remove_atoms.get_mut(atom.index()) else {
            return Err(TopologyEditError::AtomOutOfRange { atom, atom_count });
        };
        *slot = true;
        Ok(())
    }

    pub fn remove_bond(&mut self, bond: BondId) -> Result<(), TopologyEditError> {
        let bond_count = self.working.bonds.len();
        let Some(slot) = self.remove_bonds.get_mut(bond.index()) else {
            return Err(TopologyEditError::BondOutOfRange { bond, bond_count });
        };
        *slot = true;
        Ok(())
    }

    pub fn abort(self) {
        // RDKit✔️✔️: dp_delAtoms.reset();
        // RDKit✔️✔️: dp_delBonds.reset();
    }

    pub fn finish(mut self) -> Result<(TopologyBlock, TopologyMapping), TopologyEditError> {
        // RDKit✔️✔️: batchRemoveBonds();
        // RDKit✔️✔️: batchRemoveAtoms();
        for (index, bond) in self.working.bonds.iter().enumerate() {
            if self.remove_atoms[bond.begin().index()] || self.remove_atoms[bond.end().index()] {
                self.remove_bonds[index] = true;
            }
        }
        let old_atom_count = self.source.atoms.len();
        let old_bond_count = self.source.bonds.len();
        let mut atom_old_to_new = vec![None; old_atom_count];
        let mut atom_new_to_old = Vec::new();
        let mut all_atom_to_new = vec![None; self.working.atoms.len()];
        let mut next_atom_index = 0usize;
        for atom in &self.working.atoms {
            let old_index = atom.id().index();
            if self.remove_atoms[old_index] {
                continue;
            }
            let new_id = AtomId::new(next_atom_index);
            all_atom_to_new[old_index] = Some(new_id);
            if old_index < old_atom_count {
                atom_old_to_new[old_index] = Some(new_id);
                atom_new_to_old.push(Some(atom.id()));
            } else {
                atom_new_to_old.push(None);
            }
            next_atom_index += 1;
        }
        let mut atoms = Vec::with_capacity(next_atom_index);
        for atom in &self.working.atoms {
            let old_index = atom.id().index();
            let Some(new_id) = all_atom_to_new[old_index] else {
                continue;
            };
            let mut remapped = atom.clone().with_id(new_id);
            remapped
                .remap_template_attachment_order(&all_atom_to_new)
                .map_err(|source| TopologyEditError::TemplateAttachmentRemap {
                    carrier: atom.id(),
                    source,
                })?;
            atoms.push(remapped);
        }
        let mut bond_old_to_new = vec![None; old_bond_count];
        let mut bond_new_to_old = Vec::new();
        let mut bonds = Vec::new();
        for bond in &self.working.bonds {
            let old_index = bond.id().index();
            if self.remove_bonds[old_index] {
                continue;
            }
            let Some(begin) = all_atom_to_new[bond.begin().index()] else {
                continue;
            };
            let Some(end) = all_atom_to_new[bond.end().index()] else {
                continue;
            };
            let mut remapped_bond = bond.clone();
            let lost_stereo_bond = bond.stereo_atoms().is_some_and(|[left, right]| {
                self.working
                    .bonds
                    .iter()
                    .enumerate()
                    .any(|(candidate_index, candidate)| {
                        self.remove_bonds[candidate_index]
                            && (((candidate.begin() == bond.begin() && candidate.end() == left)
                                || (candidate.end() == bond.begin() && candidate.begin() == left))
                                || ((candidate.begin() == bond.end() && candidate.end() == right)
                                    || (candidate.end() == bond.end()
                                        && candidate.begin() == right)))
                    })
            });
            let stereo_atoms = (!lost_stereo_bond)
                .then(|| {
                    remapped_bond.stereo_atoms().and_then(|[left, right]| {
                        Some([
                            all_atom_to_new.get(left.index()).and_then(|x| *x)?,
                            all_atom_to_new.get(right.index()).and_then(|x| *x)?,
                        ])
                    })
                })
                .flatten();
            if stereo_atoms.is_none()
                && matches!(remapped_bond.stereo(), BondStereo::Cis | BondStereo::Trans)
            {
                // RDKit✔️✔️: if (obnd->getStereo() == Bond::BondStereo::STEREOCIS ||
                // RDKit✔️✔️:     obnd->getStereo() == Bond::BondStereo::STEREOTRANS) {
                // RDKit✔️✔️:   obnd->setStereo(Bond::BondStereo::STEREONONE);
                // RDKit✔️✔️: }
                remapped_bond
                    .set_stereo(BondStereo::None)
                    .map_err(TopologyEditError::InvalidBond)?;
            }
            let new_id = BondId::new(bonds.len());
            if old_index < old_bond_count {
                bond_old_to_new[old_index] = Some(new_id);
                bond_new_to_old.push(Some(bond.id()));
            } else {
                bond_new_to_old.push(None);
            }
            bonds.push(remapped_bond.remapped(new_id, begin, end, stereo_atoms));
        }
        let mut survives: Vec<_> = self
            .working
            .substance_groups
            .iter()
            .map(|sg| sg.can_remap_without_parent(&all_atom_to_new, &bond_old_to_new))
            .collect();
        loop {
            let mut changed = false;
            for idx in 0..self.working.substance_groups.len() {
                if survives[idx]
                    && self.working.substance_groups[idx]
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
        let mut sgroup_map = vec![None; self.working.substance_groups.len()];
        let mut next_sgroup_index = 0usize;
        for (idx, keep) in survives.iter().copied().enumerate() {
            if keep {
                sgroup_map[idx] = Some(crate::SubstanceGroupId::new(next_sgroup_index));
                next_sgroup_index += 1;
            }
        }
        let substance_groups = self
            .working
            .substance_groups
            .iter()
            .enumerate()
            .filter_map(|(idx, sg)| {
                sgroup_map[idx]
                    .and_then(|id| sg.remapped(id, &all_atom_to_new, &bond_old_to_new, &sgroup_map))
            })
            .collect();
        let stereo_groups = self
            .working
            .stereo_groups
            .iter()
            .filter_map(|group| {
                let mut group = group.clone();
                for (index, removed) in self.remove_atoms.iter().copied().enumerate() {
                    while removed && group.atoms().contains(&AtomId::new(index)) {
                        group.remove_atom(AtomId::new(index));
                    }
                }
                for (index, removed) in self.remove_bonds.iter().copied().enumerate() {
                    while removed && group.bonds().contains(&BondId::new(index)) {
                        group.remove_bond(BondId::new(index));
                    }
                }
                (!group.is_empty())
                    .then(|| group.remapped(&all_atom_to_new, &bond_old_to_new))
                    .flatten()
            })
            .collect();
        let mapping = TopologyMapping {
            atoms: AtomMapping {
                old_to_new: atom_old_to_new,
                new_to_old: atom_new_to_old,
            },
            bonds: BondMapping {
                old_to_new: bond_old_to_new,
                new_to_old: bond_new_to_old,
            },
        };
        mapping
            .validate_for_counts(old_atom_count, atoms.len(), old_bond_count, bonds.len())
            .map_err(TopologyEditError::InvalidMapping)?;
        let topology = TopologyBlock::try_from_parts(atoms, bonds, substance_groups, stereo_groups)
            .map_err(TopologyEditError::InvalidResult)?;
        Ok((topology, mapping))
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

    #[test]
    fn topology_mapping_validates_bidirectional_atom_and_bond_rows() {
        let mapping = TopologyMapping::with_appended(2, 1, 1, 1);
        assert!(mapping.validate_for_counts(2, 3, 1, 2).is_ok());

        let invalid = TopologyMapping {
            atoms: AtomMapping {
                old_to_new: vec![Some(AtomId::new(1)), Some(AtomId::new(1))],
                new_to_old: vec![None, Some(AtomId::new(0))],
            },
            bonds: BondMapping {
                old_to_new: vec![Some(BondId::new(0))],
                new_to_old: vec![Some(BondId::new(0))],
            },
        };
        assert!(matches!(
            invalid.validate_for_counts(2, 2, 1, 1),
            Err(MappingValidationError::InverseMismatch {
                entity: "atom",
                direction: "old-to-new",
                row: 1,
                mapped: 1,
            })
        ));
    }

    #[test]
    fn topology_mapping_rejects_out_of_range_targets() {
        let invalid = TopologyMapping {
            atoms: AtomMapping {
                old_to_new: vec![Some(AtomId::new(2))],
                new_to_old: vec![Some(AtomId::new(0))],
            },
            bonds: BondMapping {
                old_to_new: Vec::new(),
                new_to_old: Vec::new(),
            },
        };
        assert!(matches!(
            invalid.validate_for_counts(1, 1, 0, 0),
            Err(MappingValidationError::OutOfRange {
                entity: "atom",
                direction: "old-to-new",
                row: 0,
                mapped: 2,
                target_count: 1,
            })
        ));
    }
}
