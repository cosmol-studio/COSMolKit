//! Shared topology value container.
//!
//! This type stores concrete graph rows and local adjacency/annotation data.
//! Topology remapping policies and live-molecule installation remain owned by
//! the runtime crate.

use std::collections::HashSet;

use crate::sgroup::{SubstanceGroupValidationError, validate_substance_groups};
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

impl From<SubstanceGroupValidationError> for TopologyValidationError {
    fn from(error: SubstanceGroupValidationError) -> Self {
        match error {
            SubstanceGroupValidationError::IdMismatch { position, id } => {
                Self::SubstanceGroupIdMismatch { position, id }
            }
            SubstanceGroupValidationError::AtomOutOfRange {
                sgroup,
                atom,
                atom_count,
            } => Self::SubstanceGroupAtomOutOfRange {
                sgroup,
                atom,
                atom_count,
            },
            SubstanceGroupValidationError::BondOutOfRange {
                sgroup,
                bond,
                bond_count,
            } => Self::SubstanceGroupBondOutOfRange {
                sgroup,
                bond,
                bond_count,
            },
            SubstanceGroupValidationError::ParentOutOfRange { sgroup, parent } => {
                Self::SubstanceGroupParentOutOfRange { sgroup, parent }
            }
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, thiserror::Error)]
pub enum BondEndPointsParseErrorKind {
    #[error("no decimal prefix was found")]
    InvalidArgument,
    #[error("decimal value is outside the binary64 range")]
    OutOfRange,
    #[error("truncated decimal value is outside the u32 range")]
    UnrepresentableEndpoint,
}

#[derive(Debug, Clone, PartialEq, Eq, thiserror::Error)]
pub enum TopologyEditError {
    #[error("invalid source topology: {0}")]
    InvalidSource(TopologyValidationError),
    #[error("atom {atom} is out of range for {atom_count} atoms")]
    AtomOutOfRange { atom: AtomId, atom_count: usize },
    #[error("bond {bond} is out of range for {bond_count} bonds")]
    BondOutOfRange { bond: BondId, bond_count: usize },
    #[error("ENDPTS token {token_index} {token:?} on bond {bond} failed: {kind}")]
    BondEndPointsParse {
        bond: BondId,
        token_index: usize,
        token: String,
        kind: BondEndPointsParseErrorKind,
    },
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
        validate_substance_groups(&self.substance_groups, atom_count, bond_count)?;
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

fn parse_endpts_decimal_prefix(token: &str) -> Result<(f64, usize), BondEndPointsParseErrorKind> {
    // BEGIN RDKIT CPP FUNCTION RWMol::batchRemoveAtoms (ENDPTS conversions)
    // RDKit❗✔️: unsigned int num_ats = std::stod(*tokens.begin());
    // RDKit❗✔️: std::transform(beg, tokens.end(), std::back_inserter(oats),
    // RDKit❗✔️:                [](const std::string &a) { return std::stod(a); });
    // END RDKIT CPP FUNCTION
    // F13-NORMAL authorizes an atof-like decimal-prefix scope: ordinary finite
    // decimal values and reasonable scientific notation use Rust's binary64
    // conversion after this ASCII scan. This does not claim equivalence for
    // every std::stod/libc spelling or rare rounding boundary.
    // A scanned normal token is O(n), allocation-free; Rust parses the same
    // slice once more, so this keeps linear complexity with a second scan.
    let bytes = token.as_bytes();
    let mut index = 0;
    while index < bytes.len()
        && bytes[index] != 0
        && matches!(bytes[index], b' ' | b'\t' | b'\n' | 0x0b | 0x0c | b'\r')
    {
        index += 1;
    }
    let number_start = index;
    if index < bytes.len() && matches!(bytes[index], b'+' | b'-') {
        index += 1;
    }

    let mut digits = 0;
    while index < bytes.len() && bytes[index] != 0 && bytes[index].is_ascii_digit() {
        digits += 1;
        index += 1;
    }
    if index < bytes.len() && bytes[index] == b'.' {
        index += 1;
        while index < bytes.len() && bytes[index] != 0 && bytes[index].is_ascii_digit() {
            digits += 1;
            index += 1;
        }
    }
    if digits == 0 {
        return Err(BondEndPointsParseErrorKind::InvalidArgument);
    }

    if index < bytes.len() && matches!(bytes[index], b'e' | b'E') {
        let exponent_start = index;
        index += 1;
        if index < bytes.len() && matches!(bytes[index], b'+' | b'-') {
            index += 1;
        }
        let exponent_digits_start = index;
        while index < bytes.len() && bytes[index] != 0 && bytes[index].is_ascii_digit() {
            index += 1;
        }
        if index == exponent_digits_start {
            index = exponent_start;
        }
    }

    // Leading C whitespace contributes to std::stod's consumed index, but
    // Rust's FromStr grammar receives only the numeric slice.
    let parsed = token[number_start..index]
        .parse::<f64>()
        .map_err(|_| BondEndPointsParseErrorKind::OutOfRange)?;
    if !parsed.is_finite() {
        return Err(BondEndPointsParseErrorKind::OutOfRange);
    }
    Ok((parsed, index))
}

fn endpts_value_to_u32(value: f64) -> Result<u32, BondEndPointsParseErrorKind> {
    // This guard intentionally defines the approved safe behavior for values
    // whose C++ float-to-unsigned conversion is undefined. It is separate from
    // std::stod InvalidArgument/OutOfRange: truncate first, allow -0 and
    // fractions in (-1, 0), and cast only after the finite range check.
    if !value.is_finite() {
        return Err(BondEndPointsParseErrorKind::UnrepresentableEndpoint);
    }
    let truncated = value.trunc();
    if !(0.0..=u32::MAX as f64).contains(&truncated) {
        return Err(BondEndPointsParseErrorKind::UnrepresentableEndpoint);
    }
    Ok(truncated as u32)
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
        // RDKit❗✔️: // remove bonds attached to the atom
        // RDKit❗✔️: // In batch mode this will schedule bond removal
        // RDKit❗✔️: std::vector<std::pair<unsigned int, unsigned int>> nbrs;
        // RDKit❗✔️: ADJ_ITER b1, b2;
        // RDKit❗✔️: boost::tie(b1, b2) = getAtomNeighbors(atom);
        // RDKit❗✔️: while (b1 != b2) {
        // RDKit❗✔️:   nbrs.emplace_back(atom->getIdx(), rdcast<unsigned int>(*b1));
        // RDKit❗✔️:   ++b1;
        // RDKit❗✔️: }
        // RDKit❗✔️: for (auto &nbr : nbrs) {
        // RDKit❗✔️:   removeBond(nbr.first, nbr.second);
        // RDKit❗✔️: }
        for (index, bond) in self.working.bonds.iter().enumerate() {
            if self.remove_atoms[bond.begin().index()] || self.remove_atoms[bond.end().index()] {
                self.remove_bonds[index] = true;
            }
        }

        // RDKit✔️🔝: auto beginAtm = bnd->getBeginAtom();
        // RDKit✔️🔝: auto endAtm = bnd->getEndAtom();
        // RDKit✔️🔝: std::vector<std::vector<Atom *>> bond_atoms = {{beginAtm, endAtm},
        // RDKit✔️🔝:                                                {endAtm, beginAtm}};
        // RDKit✔️🔝: for (const auto &atoms : bond_atoms) {
        // RDKit✔️🔝:   for (auto obnd : atomBonds(atoms[0])) {
        // RDKit✔️🔝:     if (obnd == bnd) {
        // RDKit✔️🔝:       continue;
        // RDKit✔️🔝:     }
        // RDKit✔️🔝:     if (std::find(obnd->getStereoAtoms().begin(),
        // RDKit✔️🔝:                   obnd->getStereoAtoms().end(),
        // RDKit✔️🔝:                   atoms[1]->getIdx()) != obnd->getStereoAtoms().end()) {
        // RDKit✔️🔝:       if (obnd->getStereo() == Bond::BondStereo::STEREOCIS ||
        // RDKit✔️🔝:           obnd->getStereo() == Bond::BondStereo::STEREOTRANS) {
        // RDKit✔️🔝:         obnd->setStereo(Bond::BondStereo::STEREONONE);
        // RDKit✔️🔝:       }
        // RDKit✔️🔝:       obnd->getStereoAtoms().clear();
        // RDKit✔️🔝:     }
        // RDKit✔️🔝:   }
        // RDKit✔️🔝: }
        // For D scheduled bonds and S retained stereo bonds, this uses
        // expected O(D + S) index work and O(D) temporary entries. The source
        // rechecks incident neighbors for each sequential deletion; clearing
        // a stereo vector is idempotent, so indexed membership preserves the
        // final state while avoiding repeated incident scans.
        let mut removed_bond_endpoints = HashSet::new();
        for (removed, bond) in self.remove_bonds.iter().copied().zip(&self.working.bonds) {
            if !removed {
                continue;
            }
            let begin = bond.begin().index();
            let end = bond.end().index();
            removed_bond_endpoints.insert(if begin < end {
                (begin, end)
            } else {
                (end, begin)
            });
        }
        let removed_bond_connects = |left: AtomId, right: AtomId| {
            let left = left.index();
            let right = right.index();
            let pair = if left < right {
                (left, right)
            } else {
                (right, left)
            };
            removed_bond_endpoints.contains(&pair)
        };

        // BEGIN RDKIT CPP FUNCTION RWMol::batchRemoveAtoms ENDPTS updates
        // RDKit❗✔️:         if ('(' == sprop.front() && ')' == sprop.back()) {
        // RDKit❗✔️:           sprop = sprop.substr(1, sprop.length() - 2);
        // RDKit❗✔️:           boost::char_separator<char> sep(" ");
        // RDKit❗✔️:           boost::tokenizer<boost::char_separator<char>> tokens(sprop, sep);
        // RDKit❗✔️:           unsigned int num_ats = std::stod(*tokens.begin());
        // RDKit❗✔️:           std::vector<unsigned int> oats;
        // RDKit❗✔️:           auto beg = tokens.begin();
        // RDKit❗✔️:           ++beg;
        // RDKit❗✔️:           std::transform(beg, tokens.end(), std::back_inserter(oats),
        // RDKit❗✔️:                          [](const std::string &a) { return std::stod(a); });
        // RDKit❗✔️:           auto idx_pos = std::find(oats.begin(), oats.end(), idx + 1);
        // RDKit❗✔️:           if (idx_pos != oats.end()) {
        // RDKit❗✔️:             oats.erase(idx_pos);
        // RDKit❗✔️:             --num_ats;
        // RDKit❗✔️:           }
        // RDKit❗✔️:           if (!num_ats) {
        // RDKit❗✔️:             bond->clearProp(RDKit::common_properties::_MolFileBondEndPts);
        // RDKit❗✔️:             bond->clearProp(common_properties::_MolFileBondAttach);
        // RDKit❗✔️:           } else {
        // RDKit❗✔️:             sprop = "(" + std::to_string(num_ats) + " ";
        // RDKit❗✔️:             for (auto &i : oats) {
        // RDKit❗✔️:               if (i > idx + 1) {
        // RDKit❗✔️:                 --i;
        // RDKit❗✔️:               }
        // RDKit❗✔️:               sprop += std::to_string(i) + " ";
        // RDKit❗✔️:             }
        // RDKit❗✔️:             sprop[sprop.length() - 1] = ')';
        // RDKit❗✔️:             bond->setProp(RDKit::common_properties::_MolFileBondEndPts, sprop);
        // RDKit❗✔️:           }
        // RDKit❗✔️:         }
        // END RDKIT CPP FUNCTION RWMol::batchRemoveAtoms ENDPTS updates
        // Process deletions in the pinned descending order, and for each one
        // visit only bonds that survive the complete batch. Literal-space
        // splitting drops empty fields like the source tokenizer; count is
        // converted before endpoint values, then only the first matching ID
        // is erased. A missing count has no source-defined iterator result;
        // this authorized boundary returns typed InvalidArgument instead.
        // The source find/erase and endpoint remap are linear in the list;
        // this keeps the same per-deletion scan and allocates one value vector
        // and serialized property only for each parenthesized surviving bond.
        for removed_atom_index in (0..self.remove_atoms.len()).rev() {
            if !self.remove_atoms[removed_atom_index] {
                continue;
            }
            for (bond_index, bond) in self.working.bonds.iter_mut().enumerate() {
                if self.remove_bonds[bond_index] {
                    continue;
                }
                let bond_id = bond.id();
                let Some(source_property) = bond.prop("_MolFileBondEndPts") else {
                    continue;
                };
                let source_property = match source_property {
                    crate::PropertyValue::String(value) => value.clone(),
                    // RDKit✔️✔️: rdvalue_tostring(i.val, res);
                    // RDKit✔️✔️: if ('(' == sprop.front() && ')' == sprop.back()) {
                    // Dict::getValIfPresent(string&) projects Int/Double/Bool
                    // through lexical_cast. None of those scalar spellings
                    // begins with '(' (including signed zero, Inf and NaN),
                    // so the source skips this branch. Avoid that unobserved
                    // allocation without moving formatting into the model.
                    crate::PropertyValue::Int(_)
                    | crate::PropertyValue::Double(_)
                    | crate::PropertyValue::Bool(_) => continue,
                };
                let Some(contents) = source_property
                    .strip_prefix('(')
                    .and_then(|value| value.strip_suffix(')'))
                else {
                    continue;
                };

                let mut tokens = contents.split(' ').filter(|token| !token.is_empty());
                let Some(count_token) = tokens.next() else {
                    return Err(TopologyEditError::BondEndPointsParse {
                        bond: bond_id,
                        token_index: 0,
                        token: String::new(),
                        kind: BondEndPointsParseErrorKind::InvalidArgument,
                    });
                };
                let parse_token =
                    |token: &str, token_index: usize| -> Result<u32, TopologyEditError> {
                        let (value, _) = parse_endpts_decimal_prefix(token).map_err(|kind| {
                            TopologyEditError::BondEndPointsParse {
                                bond: bond_id,
                                token_index,
                                token: token.to_owned(),
                                kind,
                            }
                        })?;
                        endpts_value_to_u32(value).map_err(|kind| {
                            TopologyEditError::BondEndPointsParse {
                                bond: bond_id,
                                token_index,
                                token: token.to_owned(),
                                kind,
                            }
                        })
                    };
                let mut endpoint_count = parse_token(count_token, 0)?;
                let mut endpoints = Vec::new();
                for (endpoint_index, token) in tokens.enumerate() {
                    endpoints.push(parse_token(token, endpoint_index + 1)?);
                }

                if let Ok(removed_endpoint) = u32::try_from(removed_atom_index + 1) {
                    if let Some(position) = endpoints
                        .iter()
                        .position(|&endpoint| endpoint == removed_endpoint)
                    {
                        endpoints.remove(position);
                        // The source count is unsigned; preserve its defined
                        // wrap if a malformed count of zero still matches.
                        endpoint_count = endpoint_count.wrapping_sub(1);
                    }
                    for endpoint in &mut endpoints {
                        if *endpoint > removed_endpoint {
                            *endpoint -= 1;
                        }
                    }
                }

                if endpoint_count == 0 {
                    bond.clear_prop("_MolFileBondEndPts");
                    bond.clear_prop("_MolFileBondAttach");
                } else {
                    let mut updated = format!("({endpoint_count} ");
                    for endpoint in endpoints {
                        updated.push_str(&endpoint.to_string());
                        updated.push(' ');
                    }
                    updated.pop();
                    updated.push(')');
                    bond.set_prop("_MolFileBondEndPts", updated)
                        .map_err(TopologyEditError::InvalidBond)?;
                }
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
                removed_bond_connects(bond.begin(), left)
                    || removed_bond_connects(bond.begin(), right)
                    || removed_bond_connects(bond.end(), left)
                    || removed_bond_connects(bond.end(), right)
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
        // BEGIN RDKIT CPP FUNCTION RWMol::commitBatchEdit group-removal order
        // RDKit❗✔️: batchRemoveBonds();
        // RDKit❗✔️: batchRemoveAtoms();
        // END RDKIT CPP FUNCTION
        // BEGIN RDKIT CPP FUNCTION RWMol::batchRemoveBonds group removals
        // RDKit❗✔️:     removeSubstanceGroupsReferencingBond(*this, idx);
        // RDKit❗✔️:     removeBondFromGroups(bnd, d_stereo_groups);
        // END RDKIT CPP FUNCTION
        // BEGIN RDKIT CPP FUNCTION RWMol::batchRemoveAtoms group removals
        // RDKit❗✔️:     removeSubstanceGroupsReferencingAtom(*this, idx);
        // RDKit❗✔️:     removeAtomFromGroups(atom, d_stereo_groups);
        // END RDKIT CPP FUNCTION
        // BEGIN RDKIT CPP FUNCTION removeSubstanceGroupsReferencing helpers
        // RDKit❗✔️: bool removedParentInHierarchy(
        // RDKit❗✔️:     unsigned int idx, const std::vector<SubstanceGroup> &sgs,
        // RDKit❗✔️:     const boost::dynamic_bitset<> &toRemove,
        // RDKit❗✔️:     const std::map<unsigned int, unsigned int> &indexLookup) {
        // RDKit❗✔️:   PRECONDITION(idx < sgs.size(), "cannot find SubstanceGroup");
        // RDKit❗✔️:   if (toRemove[idx]) {
        // RDKit❗✔️:     return true;
        // RDKit❗✔️:   }
        // RDKit❗✔️:   unsigned int parent;
        // RDKit❗✔️:   if (sgs[idx].getPropIfPresent("PARENT", parent)) {
        // RDKit❗✔️:     auto piter = indexLookup.find(parent);
        // RDKit❗✔️:     if (piter != indexLookup.end()) {
        // RDKit❗✔️:       return removedParentInHierarchy(piter->second, sgs, toRemove,
        // RDKit❗✔️:                                       indexLookup);
        // RDKit❗✔️:     }
        // RDKit❗✔️:   }
        // RDKit❗✔️:   return false;
        // RDKit❗✔️: }
        // RDKit❗✔️: template <bool INCLUDES_METHOD(SubstanceGroup &, unsigned int),
        // RDKit❗✔️:           void ADJUST_METHOD(SubstanceGroup &, unsigned int)>
        // RDKit❗✔️: void removeSubstanceGroupsReferencing(RWMol &mol, unsigned int idx) {
        // RDKit❗✔️:   auto &sgs = getSubstanceGroups(mol);
        // RDKit❗✔️:   if (!sgs.empty()) {
        // RDKit❗✔️:     // first collect the ones that should be removed
        // RDKit❗✔️:     boost::dynamic_bitset<> toRemove(sgs.size());
        // RDKit❗✔️:     unsigned int nRemoved = 0;
        // RDKit❗✔️:     bool parentsPresent = false;
        // RDKit❗✔️:     for (unsigned int i = 0; i < sgs.size(); ++i) {
        // RDKit❗✔️:       if (!parentsPresent && sgs[i].hasProp("PARENT")) {
        // RDKit❗✔️:         parentsPresent = true;
        // RDKit❗✔️:       }
        // RDKit❗✔️:       if (INCLUDES_METHOD(sgs[i], idx)) {
        // RDKit❗✔️:         toRemove.set(i);
        // RDKit❗✔️:         ++nRemoved;
        // RDKit❗✔️:       }
        // RDKit❗✔️:     }
        // RDKit❗✔️:     // if we're going to be removing anything and there are PARENTS present,
        // RDKit❗✔️:     // we need to build a lookup map between index->position in original array
        // RDKit❗✔️:     std::map<unsigned int, unsigned int> indexLookup;
        // RDKit❗✔️:     if (parentsPresent && nRemoved) {
        // RDKit❗✔️:       for (unsigned int i = 0; i < sgs.size(); ++i) {
        // RDKit❗✔️:         unsigned int index;
        // RDKit❗✔️:         if (sgs[i].getPropIfPresent("index", index)) {
        // RDKit❗✔️:           indexLookup[index] = i;
        // RDKit❗✔️:         }
        // RDKit❗✔️:       }
        // RDKit❗✔️:     }
        // RDKit❗✔️:     // now go through and keep everything that shouldn't be removed
        // RDKit❗✔️:     // and who doesn't have a PARENT that should be removed in their hierarchy
        // RDKit❗✔️:     std::vector<SubstanceGroup> newsgs;
        // RDKit❗✔️:     newsgs.reserve(sgs.size() - nRemoved);
        // RDKit❗✔️:     unsigned int i = 0;
        // RDKit❗✔️:     for (auto &&sg : sgs) {
        // RDKit❗✔️:       if (!toRemove[i]) {
        // RDKit❗✔️:         // we might be keeping it. Check the parent
        // RDKit❗✔️:         if (!parentsPresent || !sg.hasProp("PARENT")) {
        // RDKit❗✔️:           ADJUST_METHOD(sg, idx);
        // RDKit❗✔️:           newsgs.push_back(std::move(sg));
        // RDKit❗✔️:         } else if (parentsPresent) {
        // RDKit❗✔️:           unsigned int parent;
        // RDKit❗✔️:           // has our parent been removed?
        // RDKit❗✔️:           if (sg.getPropIfPresent("PARENT", parent)) {
        // RDKit❗✔️:             auto piter = indexLookup.find(parent);
        // RDKit❗✔️:             bool keepIt = false;
        // RDKit❗✔️:             if (piter == indexLookup.end()) {
        // RDKit❗✔️:               // our parent isn't around, so it isn't being removed
        // RDKit❗✔️:               // note: this is an odd case and probably shouldn't happen, but
        // RDKit❗✔️:               // this isn't the place to enforce that
        // RDKit❗✔️:               keepIt = true;
        // RDKit❗✔️:             } else if (!toRemove[piter->second]) {
        // RDKit❗✔️:               // our parent isn't being removed, recursively check up through
        // RDKit❗✔️:               // parents to see if we find any that are being removed:
        // RDKit❗✔️:               if (!removedParentInHierarchy(piter->second, sgs, toRemove,
        // RDKit❗✔️:                                             indexLookup)) {
        // RDKit❗✔️:                 keepIt = true;
        // RDKit❗✔️:               }
        // RDKit❗✔️:             }
        // RDKit❗✔️:             if (keepIt) {
        // RDKit❗✔️:               ADJUST_METHOD(sg, idx);
        // RDKit❗✔️:               newsgs.push_back(std::move(sg));
        // RDKit❗✔️:             }
        // RDKit❗✔️:           }
        // RDKit❗✔️:         }
        // RDKit❗✔️:       }
        // RDKit❗✔️:       ++i;
        // RDKit❗✔️:     }
        // RDKit❗✔️:     sgs = std::move(newsgs);
        // RDKit❗✔️:   }
        // RDKit❗✔️: }
        // RDKit❗✔️: void removeSubstanceGroupsReferencingAtom(RWMol &mol, unsigned int idx) {
        // RDKit❗✔️:   // Delete substance groups containing this atom. It could be that it's ok to
        // RDKit❗✔️:   // keep it, but we just don't know
        // RDKit❗✔️:   removeSubstanceGroupsReferencing<includesAtom, removedAtom>(mol, idx);
        // RDKit❗✔️: }
        // RDKit❗✔️: void removeSubstanceGroupsReferencingBond(RWMol &mol, unsigned int idx) {
        // RDKit❗✔️:   // Delete substance groups containing this bond. It could be that it's ok to
        // RDKit❗✔️:   // keep it, but we just don't know
        // RDKit❗✔️:   removeSubstanceGroupsReferencing<includesBond, removedBond>(mol, idx);
        // RDKit❗✔️: }
        // END RDKIT CPP FUNCTION
        // BEGIN RDKIT CPP FUNCTION StereoGroup::removeBondFromGroups
        // RDKit❗✔️: void removeBondFromGroups(const Bond *bond, std::vector<StereoGroup> &groups) {
        // RDKit❗✔️:   auto findBond = [bond](StereoGroup &group) {
        // RDKit❗✔️:     return std::find(group.getBonds().begin(), group.getBonds().end(), bond);
        // RDKit❗✔️:   };
        // RDKit❗✔️:   for (auto &group : groups) {
        // RDKit❗✔️:     auto bondPos = findBond(group);
        // RDKit❗✔️:     if (bondPos != group.d_bonds.end()) {
        // RDKit❗✔️:       group.d_bonds.erase(bondPos);
        // RDKit❗✔️:     }
        // RDKit❗✔️:   }
        // RDKit❗✔️:   // now remove any empty groups:
        // RDKit❗✔️:   groups.erase(std::remove_if(groups.begin(), groups.end(),
        // RDKit❗✔️:                               [](const auto &gp) {
        // RDKit❗✔️:                                 return gp.getAtoms().empty() &&
        // RDKit❗✔️:                                        gp.getBonds().empty();
        // RDKit❗✔️:                               }),
        // RDKit❗✔️:                groups.end());
        // RDKit❗✔️: }
        // END RDKIT CPP FUNCTION
        // BEGIN RDKIT CPP FUNCTION StereoGroup::removeAtomFromGroups
        // RDKit❗✔️: void removeAtomFromGroups(const Atom *atom, std::vector<StereoGroup> &groups) {
        // RDKit❗✔️:   auto findAtom = [atom](StereoGroup &group) {
        // RDKit❗✔️:     return std::find(group.getAtoms().begin(), group.getAtoms().end(), atom);
        // RDKit❗✔️:   };
        // RDKit❗✔️:   for (auto &group : groups) {
        // RDKit❗✔️:     auto atomPos = findAtom(group);
        // RDKit❗✔️:     if (atomPos != group.d_atoms.end()) {
        // RDKit❗✔️:       group.d_atoms.erase(atomPos);
        // RDKit❗✔️:     }
        // RDKit❗✔️:   }
        // RDKit❗✔️:   // now remove any empty groups:
        // RDKit❗✔️:   groups.erase(std::remove_if(groups.begin(), groups.end(),
        // RDKit❗✔️:                               [](const auto &gp) {
        // RDKit❗✔️:                                 return gp.getAtoms().empty() &&
        // RDKit❗✔️:                                        gp.getBonds().empty();
        // RDKit❗✔️:                               }),
        // RDKit❗✔️:                groups.end());
        // RDKit❗✔️: }
        // END RDKIT CPP FUNCTION
        // The source helpers scan in their stated batch order and preserve group
        // and surviving-member order. The detached maps below combine the same
        // row shifts, while this source-shaped member pass removes one matching
        // entry for each removed ID before compacting the groups.
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
            .filter_map(|source_group| {
                let mut group = source_group.clone();
                for (index, removed) in self.remove_bonds.iter().copied().enumerate().rev() {
                    if removed {
                        group.remove_bond(BondId::new(index));
                    }
                }
                for (index, removed) in self.remove_atoms.iter().copied().enumerate().rev() {
                    if removed {
                        group.remove_atom(AtomId::new(index));
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
    use crate::{
        AtomSpec, BondOrder, BondSpec, SGroupAttachPoint, SGroupBondRole, SGroupCState,
        StereoGroupKind, SubstanceGroupId, SubstanceGroupKind,
    };

    fn atom(id: usize) -> Atom {
        Atom::from_spec(AtomId::new(id), AtomSpec::new(crate::Element::C))
    }

    fn f12_bond(
        id: usize,
        begin: usize,
        end: usize,
        order: BondOrder,
        stereo: BondStereo,
        stereo_atoms: Option<[usize; 2]>,
    ) -> Bond {
        let mut spec =
            BondSpec::new(AtomId::new(begin), AtomId::new(end), order).with_stereo(stereo);
        if let Some([left, right]) = stereo_atoms {
            spec = spec.with_stereo_atoms(AtomId::new(left), AtomId::new(right));
        }
        Bond::from_spec(BondId::new(id), spec)
    }

    fn f12_topology(atom_count: usize, bonds: Vec<Bond>) -> TopologyBlock {
        TopologyBlock::try_from_parts(
            (0..atom_count).map(atom).collect(),
            bonds,
            Vec::new(),
            Vec::new(),
        )
        .expect("fixed F12 topology is valid")
    }

    fn f13_bond(
        id: usize,
        begin: usize,
        end: usize,
        endpts: Option<&str>,
        attach: Option<&str>,
    ) -> Bond {
        let mut spec = BondSpec::new(AtomId::new(begin), AtomId::new(end), BondOrder::Single);
        if let Some(value) = endpts {
            spec = spec
                .with_prop("_MolFileBondEndPts", value)
                .expect("fixed ENDPTS key is valid");
        }
        if let Some(value) = attach {
            spec = spec
                .with_prop("_MolFileBondAttach", value)
                .expect("fixed ATTACH key is valid");
        }
        Bond::from_spec(BondId::new(id), spec)
    }

    fn f14_topology() -> TopologyBlock {
        let atoms: Vec<_> = (0..8).map(atom).collect();
        let bonds = vec![
            f12_bond(0, 0, 1, BondOrder::Single, BondStereo::None, None),
            f12_bond(1, 2, 3, BondOrder::Single, BondStereo::None, None),
            f12_bond(2, 5, 6, BondOrder::Single, BondStereo::None, None),
            f12_bond(3, 6, 7, BondOrder::Single, BondStereo::None, None),
            f12_bond(4, 0, 7, BondOrder::Single, BondStereo::None, None),
        ];

        let mut root = SubstanceGroup::new(SubstanceGroupId::new(0), SubstanceGroupKind::Data)
            .with_rdkit_sequence_id(71)
            .with_external_id(501)
            .with_atoms(vec![AtomId::new(7), AtomId::new(0), AtomId::new(5)])
            .with_bonds(vec![BondId::new(3), BondId::new(0), BondId::new(2)])
            .with_bond_role(BondId::new(3), SGroupBondRole::Contained)
            .with_head_crossing_bonds(vec![BondId::new(3), BondId::new(2)])
            .with_crossing_bond_correspondence(vec![BondId::new(2), BondId::new(0)])
            .with_parent_atoms(vec![AtomId::new(6)])
            .with_attach_points(vec![SGroupAttachPoint {
                atom: AtomId::new(5),
                leaving_atom: Some(AtomId::new(7)),
                label: Some("AP".into()),
                order: Some(4),
            }])
            .with_cstates(vec![SGroupCState::new(BondId::new(2), [1.0, 2.0, 3.0])])
            .with_label("retained-root")
            .with_data_field("root-field");
        root.set_prop("vendor", "keep");

        let groups = vec![
            root,
            SubstanceGroup::new(SubstanceGroupId::new(1), SubstanceGroupKind::Data)
                .with_atoms(vec![AtomId::new(2)]),
            SubstanceGroup::new(SubstanceGroupId::new(2), SubstanceGroupKind::Data)
                .with_atoms(vec![AtomId::new(4)])
                .with_parent(SubstanceGroupId::new(1)),
            SubstanceGroup::new(SubstanceGroupId::new(3), SubstanceGroupKind::Data)
                .with_atoms(vec![AtomId::new(6), AtomId::new(0)])
                .with_bonds(vec![BondId::new(3), BondId::new(2)])
                .with_parent(SubstanceGroupId::new(0)),
            SubstanceGroup::new(SubstanceGroupId::new(4), SubstanceGroupKind::Data)
                .with_parent_atoms(vec![AtomId::new(2)]),
            SubstanceGroup::new(SubstanceGroupId::new(5), SubstanceGroupKind::Data)
                .with_attach_points(vec![SGroupAttachPoint {
                    atom: AtomId::new(4),
                    leaving_atom: Some(AtomId::new(2)),
                    label: Some("REMOVE".into()),
                    order: Some(1),
                }]),
            SubstanceGroup::new(SubstanceGroupId::new(6), SubstanceGroupKind::Data)
                .with_bonds(vec![BondId::new(4)]),
            SubstanceGroup::new(SubstanceGroupId::new(7), SubstanceGroupKind::Data)
                .with_atoms(vec![AtomId::new(0)])
                .with_parent(SubstanceGroupId::new(6)),
            SubstanceGroup::new(SubstanceGroupId::new(8), SubstanceGroupKind::Data)
                .with_cstates(vec![SGroupCState::new(BondId::new(4), [4.0, 5.0, 6.0])]),
        ];
        let stereo_groups = vec![
            crate::StereoGroup::new(
                StereoGroupKind::And,
                vec![AtomId::new(7), AtomId::new(0), AtomId::new(5)],
                vec![BondId::new(3), BondId::new(0), BondId::new(2)],
            )
            .with_id(41)
            .with_write_id(9),
            crate::StereoGroup::new(
                StereoGroupKind::Or,
                vec![AtomId::new(5), AtomId::new(2), AtomId::new(0)],
                vec![BondId::new(2), BondId::new(4), BondId::new(0)],
            )
            .with_id(42)
            .with_write_id(12),
            crate::StereoGroup::new(StereoGroupKind::Absolute, vec![AtomId::new(2)], vec![])
                .with_id(43)
                .with_write_id(13),
            crate::StereoGroup::new(StereoGroupKind::Or, vec![], vec![BondId::new(4)])
                .with_id(44)
                .with_write_id(14),
            crate::StereoGroup::new(
                StereoGroupKind::And,
                vec![AtomId::new(2)],
                vec![BondId::new(4)],
            )
            .with_id(45)
            .with_write_id(15),
        ];

        TopologyBlock::try_from_parts(atoms, bonds, groups, stereo_groups)
            .expect("fixed F14 topology is valid")
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
    fn query_sgroups_shared_validation_preserves_topology_error_categories() {
        let id = SubstanceGroupId::new(0);
        let invalid_atom = AtomId::new(2);
        let invalid_bond = BondId::new(1);
        let atoms = vec![atom(0), atom(1)];
        let bonds = vec![Bond::from_spec(
            BondId::new(0),
            BondSpec::new(AtomId::new(0), AtomId::new(1), BondOrder::Single),
        )];
        let validate = |group| {
            TopologyBlock::try_from_parts(atoms.clone(), bonds.clone(), vec![group], Vec::new())
        };

        assert_eq!(
            validate(SubstanceGroup::new(
                SubstanceGroupId::new(1),
                SubstanceGroupKind::Data,
            )),
            Err(TopologyValidationError::SubstanceGroupIdMismatch {
                position: 0,
                id: SubstanceGroupId::new(1),
            })
        );
        for invalid in [
            SubstanceGroup::new(id, SubstanceGroupKind::Data).with_atoms(vec![invalid_atom]),
            SubstanceGroup::new(id, SubstanceGroupKind::Data).with_parent_atoms(vec![invalid_atom]),
            SubstanceGroup::new(id, SubstanceGroupKind::Data).with_attach_points(vec![
                SGroupAttachPoint {
                    atom: invalid_atom,
                    leaving_atom: None,
                    label: None,
                    order: None,
                },
            ]),
            SubstanceGroup::new(id, SubstanceGroupKind::Data).with_attach_points(vec![
                SGroupAttachPoint {
                    atom: AtomId::new(0),
                    leaving_atom: Some(invalid_atom),
                    label: None,
                    order: None,
                },
            ]),
        ] {
            assert_eq!(
                validate(invalid),
                Err(TopologyValidationError::SubstanceGroupAtomOutOfRange {
                    sgroup: id,
                    atom: invalid_atom,
                    atom_count: 2,
                })
            );
        }
        for invalid in [
            SubstanceGroup::new(id, SubstanceGroupKind::Data).with_bonds(vec![invalid_bond]),
            SubstanceGroup::new(id, SubstanceGroupKind::Data)
                .with_cstates(vec![SGroupCState::new(invalid_bond, [1.0, 2.0, 3.0])]),
            SubstanceGroup::new(id, SubstanceGroupKind::Data)
                .with_head_crossing_bonds(vec![invalid_bond]),
            SubstanceGroup::new(id, SubstanceGroupKind::Data)
                .with_crossing_bond_correspondence(vec![invalid_bond]),
        ] {
            assert_eq!(
                validate(invalid),
                Err(TopologyValidationError::SubstanceGroupBondOutOfRange {
                    sgroup: id,
                    bond: invalid_bond,
                    bond_count: 1,
                })
            );
        }
        assert_eq!(
            validate(
                SubstanceGroup::new(id, SubstanceGroupKind::Data)
                    .with_parent(SubstanceGroupId::new(1)),
            ),
            Err(TopologyValidationError::SubstanceGroupParentOutOfRange {
                sgroup: id,
                parent: SubstanceGroupId::new(1),
            })
        );
        assert!(
            TopologyBlock::try_from_parts(
                atoms,
                bonds,
                vec![
                    SubstanceGroup::new(id, SubstanceGroupKind::Data)
                        .with_parent(SubstanceGroupId::new(1)),
                    SubstanceGroup::new(SubstanceGroupId::new(1), SubstanceGroupKind::Data)
                        .with_parent(id),
                ],
                Vec::new(),
            )
            .is_ok()
        );
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

    #[test]
    fn cf3d_frag_f12_whole_component_removal_preserves_retained_stereo_and_gapped_bond_ids() {
        let source = f12_topology(
            6,
            vec![
                f12_bond(0, 0, 1, BondOrder::Double, BondStereo::Cis, Some([2, 3])),
                f12_bond(1, 4, 5, BondOrder::Single, BondStereo::None, None),
                f12_bond(2, 0, 2, BondOrder::Single, BondStereo::None, None),
                f12_bond(3, 1, 3, BondOrder::Single, BondStereo::None, None),
            ],
        );
        let original = source.clone();
        let mut edit = source.begin_batch_edit().expect("valid source topology");
        edit.remove_atom(AtomId::new(4)).unwrap();
        edit.remove_atom(AtomId::new(5)).unwrap();

        let (result, mapping) = edit.finish().expect("component deletion succeeds");

        assert_eq!(source, original);
        assert_eq!(result.bonds.len(), 3);
        assert_eq!(
            mapping.bonds.old_to_new,
            vec![
                Some(BondId::new(0)),
                None,
                Some(BondId::new(1)),
                Some(BondId::new(2)),
            ]
        );
        assert_eq!(result.bonds[0].stereo(), BondStereo::Cis);
        assert_eq!(
            result.bonds[0].stereo_atoms(),
            Some([AtomId::new(2), AtomId::new(3)])
        );
        assert_eq!(result.bonds[1].begin(), AtomId::new(0));
        assert_eq!(result.bonds[1].end(), AtomId::new(2));
        assert!(result.validate().is_ok());
    }

    #[test]
    fn cf3d_frag_f12_removed_reference_edge_checks_all_endpoint_pairs_and_stereo_tags() {
        let stereo_values = [
            BondStereo::None,
            BondStereo::Any,
            BondStereo::Z,
            BondStereo::E,
            BondStereo::Cis,
            BondStereo::Trans,
            BondStereo::AtropCw,
            BondStereo::AtropCcw,
        ];
        // Each deleted edge is reversed in its row to check undirected
        // endpoint indexing as well as every target-end/reference pairing.
        let reference_edges = [(2, 0), (3, 0), (2, 1), (3, 1)];

        for stereo in stereo_values {
            for (left, right) in reference_edges {
                let source = f12_topology(
                    5,
                    vec![
                        f12_bond(0, 0, 1, BondOrder::Double, stereo, Some([2, 3])),
                        f12_bond(1, left, right, BondOrder::Single, BondStereo::None, None),
                    ],
                );
                let mut edit = source.begin_batch_edit().expect("valid source topology");
                edit.remove_bond(BondId::new(1)).unwrap();
                let (result, _) = edit.finish().expect("reference edge removal succeeds");

                assert_eq!(result.bonds[0].stereo_atoms(), None);
                assert_eq!(
                    result.bonds[0].stereo(),
                    if matches!(stereo, BondStereo::Cis | BondStereo::Trans) {
                        BondStereo::None
                    } else {
                        stereo
                    },
                    "source tag transition for {stereo:?} after removing {left}-{right}"
                );
                assert!(result.validate().is_ok());
            }
        }
    }

    #[test]
    fn cf3d_frag_f12_removing_stereo_component_keeps_other_component_and_compacts_ids() {
        let source = f12_topology(
            7,
            vec![
                f12_bond(0, 0, 1, BondOrder::Double, BondStereo::Trans, Some([2, 3])),
                f12_bond(1, 0, 2, BondOrder::Single, BondStereo::None, None),
                f12_bond(2, 1, 3, BondOrder::Single, BondStereo::None, None),
                f12_bond(3, 4, 5, BondOrder::Single, BondStereo::None, None),
            ],
        );
        let original = source.clone();
        let mut edit = source.begin_batch_edit().expect("valid source topology");
        for atom in 0..4 {
            edit.remove_atom(AtomId::new(atom)).unwrap();
        }

        let (result, mapping) = edit.finish().expect("stereo component deletion succeeds");

        assert_eq!(source, original);
        assert_eq!(result.atoms.len(), 3);
        assert_eq!(result.bonds.len(), 1);
        assert_eq!(
            mapping.bonds.old_to_new,
            vec![None, None, None, Some(BondId::new(0))]
        );
        assert_eq!(result.bonds[0].begin(), AtomId::new(0));
        assert_eq!(result.bonds[0].end(), AtomId::new(1));
        assert!(result.validate().is_ok());
    }

    #[test]
    fn cf3d_frag_f12_noop_removal_preserves_stereo_and_identity_mapping() {
        let source = f12_topology(
            4,
            vec![
                f12_bond(0, 0, 1, BondOrder::Double, BondStereo::E, Some([2, 3])),
                f12_bond(1, 0, 2, BondOrder::Single, BondStereo::None, None),
                f12_bond(2, 1, 3, BondOrder::Single, BondStereo::None, None),
            ],
        );
        let edit = source.begin_batch_edit().expect("valid source topology");

        let (result, mapping) = edit.finish().expect("no-op finish succeeds");

        assert_eq!(result, source);
        assert_eq!(
            mapping.bonds.old_to_new,
            vec![
                Some(BondId::new(0)),
                Some(BondId::new(1)),
                Some(BondId::new(2)),
            ]
        );
        assert_eq!(result.bonds[0].stereo(), BondStereo::E);
        assert_eq!(
            result.bonds[0].stereo_atoms(),
            Some([AtomId::new(2), AtomId::new(3)])
        );
    }

    #[test]
    fn cf3d_endpts_numeric_decimal_prefix_matches_fixed_reference() {
        // Fixed in-scope values/bits and consumed-prefix lengths from the
        // ordinary std::stod reference manifest; expectations are not derived
        // from this implementation.
        let cases = [
            ("0", 0x0000_0000_0000_0000, 1),
            ("1", 0x3ff0_0000_0000_0000, 1),
            ("+1", 0x3ff0_0000_0000_0000, 2),
            ("-0", 0x8000_0000_0000_0000, 2),
            ("-0.0", 0x8000_0000_0000_0000, 4),
            (".5", 0x3fe0_0000_0000_0000, 2),
            ("1.", 0x3ff0_0000_0000_0000, 2),
            ("1.5", 0x3ff8_0000_0000_0000, 3),
            ("-0.5", 0xbfe0_0000_0000_0000, 4),
            ("-1", 0xbff0_0000_0000_0000, 2),
            ("1e2", 0x4059_0000_0000_0000, 3),
            ("1e-2", 0x3f84_7ae1_47ae_147b, 4),
            ("-1.25e-1", 0xbfc0_0000_0000_0000, 8),
            ("2.e1", 0x4034_0000_0000_0000, 4),
            ("001.0", 0x3ff0_0000_0000_0000, 5),
            ("1.25e+2kg", 0x405f_4000_0000_0000, 7),
            ("1junk", 0x3ff0_0000_0000_0000, 1),
            ("1e", 0x3ff0_0000_0000_0000, 1),
            ("1e+", 0x3ff0_0000_0000_0000, 1),
            ("1e-", 0x3ff0_0000_0000_0000, 1),
            ("1.25e+x", 0x3ff4_0000_0000_0000, 4),
            ("1\0e2", 0x3ff0_0000_0000_0000, 1),
            ("4294967295", 0x41ef_ffff_ffe0_0000, 10),
            ("4294967295.5", 0x41ef_ffff_fff0_0000, 12),
            ("4294967295.9999995", 0x41ef_ffff_ffff_ffff, 18),
            ("4294967296", 0x41f0_0000_0000_0000, 10),
        ];

        for (token, expected_bits, expected_consumed) in cases {
            let (value, consumed) =
                parse_endpts_decimal_prefix(token).expect("fixed decimal prefix converts");
            assert_eq!(
                value.to_bits(),
                expected_bits,
                "binary64 result for {token:?}"
            );
            assert_eq!(consumed, expected_consumed, "consumed prefix for {token:?}");
        }
    }

    #[test]
    fn cf3d_endpts_numeric_c_whitespace_and_no_conversion_are_fixed() {
        for (token, expected_consumed) in [
            (" 1", 2),
            ("\t1", 2),
            ("\n1", 2),
            ("\x0b1", 2),
            ("\x0c1", 2),
            ("\r1", 2),
        ] {
            let (value, consumed) =
                parse_endpts_decimal_prefix(token).expect("C whitespace precedes a number");
            assert_eq!(value.to_bits(), 0x3ff0_0000_0000_0000, "{token:?}");
            assert_eq!(consumed, expected_consumed, "{token:?}");
        }

        for token in ["", " ", "\t\n\x0b\x0c\r", "abc", "+", "-", ".", "e1", "\0"] {
            assert_eq!(
                parse_endpts_decimal_prefix(token),
                Err(BondEndPointsParseErrorKind::InvalidArgument),
                "no-conversion classification for {token:?}"
            );
        }
        assert_eq!(
            parse_endpts_decimal_prefix("1e309"),
            Err(BondEndPointsParseErrorKind::OutOfRange)
        );
    }

    #[test]
    fn cf3d_endpts_numeric_truncation_guard_has_fixed_boundaries() {
        for (value, expected) in [
            (-0.0, 0),
            (-0.5, 0),
            (f64::from_bits(0xbfefffffffffffff), 0),
            (0.0, 0),
            (u32::MAX as f64, u32::MAX),
            (u32::MAX as f64 + 0.5, u32::MAX),
            (f64::from_bits(0x41efffffffffffff), u32::MAX),
        ] {
            assert_eq!(endpts_value_to_u32(value), Ok(expected), "value {value:?}");
        }

        for value in [
            -1.0,
            f64::from_bits(0xbff0000000000001),
            4_294_967_296.0,
            f64::NAN,
            f64::INFINITY,
            f64::NEG_INFINITY,
        ] {
            assert_eq!(
                endpts_value_to_u32(value),
                Err(BondEndPointsParseErrorKind::UnrepresentableEndpoint),
                "value {value:?}"
            );
        }
    }

    #[test]
    fn cf3d_endpts_numeric_excluded_diagnostics_keep_local_behavior_explicit() {
        // These are deliberately not source-equivalence assertions. The pinned
        // std::stod oracle consumes 0x1p+2 as hexadecimal 4, while this
        // authorized decimal-only scan deterministically consumes its "0"
        // prefix. The oracle classifies 1e-310 as OutOfRange; Rust's ordinary
        // f64 conversion accepts a finite subnormal here. Hexadecimal and
        // subnormal/underflow parity remain outside F13-NORMAL.
        let (hex_prefix, hex_consumed) = parse_endpts_decimal_prefix("0x1p+2").unwrap();
        assert_eq!(hex_prefix.to_bits(), 0x0000_0000_0000_0000);
        assert_eq!(hex_consumed, 1);

        let (rust_subnormal, subnormal_consumed) = parse_endpts_decimal_prefix("1e-310").unwrap();
        assert!(rust_subnormal.is_finite() && rust_subnormal > 0.0);
        assert_eq!(subnormal_consumed, 6);
    }

    #[test]
    fn cf3d_frag_f13_absent_ids_zero_cleanup_mismatch_and_raw_values_follow_source() {
        let source = f12_topology(
            7,
            vec![
                f13_bond(0, 0, 1, Some("(2 4 5)"), Some("ANY")),
                f13_bond(1, 5, 6, Some("(0 2 7)"), Some("ANY")),
                f13_bond(2, 0, 6, Some("raw 4 5"), Some("ANY")),
                f13_bond(3, 1, 5, Some("(2.9 2.9 4.2)"), Some("ANY")),
                f13_bond(4, 3, 4, Some("(9 3 8)"), Some("ANY")),
            ],
        );
        let original = source.clone();
        let mut edit = source.begin_batch_edit().expect("valid fixed topology");
        edit.remove_atom(AtomId::new(2)).expect("atom 2 exists");

        let (result, _) = edit.finish().expect("ordinary ENDPTS values remap");

        assert_eq!(source, original, "batch edit does not mutate its source");
        assert_eq!(
            result.bonds[0].prop("_MolFileBondEndPts"),
            Some(&crate::PropertyValue::String("(2 3 4)".to_owned()))
        );
        assert_eq!(
            result.bonds[0].prop("_MolFileBondAttach"),
            Some(&crate::PropertyValue::String("ANY".to_owned()))
        );
        assert_eq!(result.bonds[1].prop("_MolFileBondEndPts"), None);
        assert_eq!(result.bonds[1].prop("_MolFileBondAttach"), None);
        assert_eq!(
            result.bonds[2].prop("_MolFileBondEndPts"),
            Some(&crate::PropertyValue::String("raw 4 5".to_owned()))
        );
        assert_eq!(
            result.bonds[2].prop("_MolFileBondAttach"),
            Some(&crate::PropertyValue::String("ANY".to_owned()))
        );
        assert_eq!(
            result.bonds[3].prop("_MolFileBondEndPts"),
            Some(&crate::PropertyValue::String("(2 2 3)".to_owned()))
        );
        assert_eq!(
            result.bonds[3].prop("_MolFileBondAttach"),
            Some(&crate::PropertyValue::String("ANY".to_owned()))
        );
        assert_eq!(
            result.bonds[4].prop("_MolFileBondEndPts"),
            Some(&crate::PropertyValue::String("(8 7)".to_owned()))
        );
        assert_eq!(
            result.bonds[4].prop("_MolFileBondAttach"),
            Some(&crate::PropertyValue::String("ANY".to_owned()))
        );
    }

    #[test]
    fn uff_sync_endpts_typed_scalars_skip_rewrite_and_keep_exact_values() {
        let scalar_values = [
            (
                "nonparenthesized String",
                crate::PropertyValue::String("raw 4 5".to_owned()),
            ),
            ("zero Int", crate::PropertyValue::Int(0)),
            ("negative Int", crate::PropertyValue::Int(-17)),
            ("signed-zero Double", crate::PropertyValue::Double(-0.0)),
            ("finite Double", crate::PropertyValue::Double(2.75)),
            ("false Bool", crate::PropertyValue::Bool(false)),
            ("true Bool", crate::PropertyValue::Bool(true)),
        ];

        for (case, endpts_value) in scalar_values {
            let typed_bond = |id, begin, end| {
                let spec = BondSpec::new(AtomId::new(begin), AtomId::new(end), BondOrder::Single)
                    .with_prop("_MolFileBondEndPts", endpts_value.clone())
                    .expect("fixed ENDPTS property is valid")
                    .with_prop("_MolFileBondAttach", "ANY")
                    .expect("fixed ATTACH property is valid");
                Bond::from_spec(BondId::new(id), spec)
            };
            let source = f12_topology(4, vec![typed_bond(0, 0, 1), typed_bond(1, 2, 3)]);
            let original = source.clone();
            let mut edit = source.begin_batch_edit().expect("valid fixed topology");
            edit.remove_atom(AtomId::new(0))
                .expect("selected atom 0 exists");

            let (result, mapping) = edit.finish().expect("typed scalar is skipped");

            assert_eq!(source, original, "source remains unchanged for {case}");
            assert_eq!(
                mapping.atoms.old_to_new(),
                &[
                    None,
                    Some(AtomId::new(0)),
                    Some(AtomId::new(1)),
                    Some(AtomId::new(2))
                ],
                "selected and unselected atom rows for {case}"
            );
            assert_eq!(
                mapping.bonds.old_to_new(),
                &[None, Some(BondId::new(0))],
                "selected-atom bond is deleted and the other survives for {case}"
            );
            assert_eq!(result.bonds.len(), 1, "surviving bond count for {case}");
            assert_eq!(result.bonds[0].id(), BondId::new(0));
            assert_eq!(
                result.bonds[0].prop("_MolFileBondEndPts"),
                Some(&endpts_value),
                "surviving ENDPTS retains its exact type and value for {case}"
            );
            assert_eq!(
                result.bonds[0].prop("_MolFileBondAttach"),
                Some(&crate::PropertyValue::String("ANY".to_owned())),
                "scalar skip leaves ATTACH unchanged for {case}"
            );
        }

        let source = f12_topology(
            4,
            vec![
                f13_bond(0, 0, 1, Some("(2 2 3)"), Some("ANY")),
                f13_bond(1, 2, 3, Some("(2 2 3)"), Some("ANY")),
            ],
        );
        let mut edit = source.begin_batch_edit().expect("valid fixed topology");
        edit.remove_atom(AtomId::new(0))
            .expect("selected atom 0 exists");

        let (result, mapping) = edit.finish().expect("parenthesized String remaps");

        assert_eq!(mapping.bonds.old_to_new(), &[None, Some(BondId::new(0))]);
        assert_eq!(
            result.bonds[0].prop("_MolFileBondEndPts"),
            Some(&crate::PropertyValue::String("(2 1 2)".to_owned()))
        );
        assert_eq!(
            result.bonds[0].prop("_MolFileBondAttach"),
            Some(&crate::PropertyValue::String("ANY".to_owned()))
        );
    }

    #[test]
    fn cf3d_frag_f13_descending_deletions_remove_only_first_duplicate() {
        let source = f12_topology(7, vec![f13_bond(0, 0, 1, Some("(4 3 3 5 6)"), Some("ANY"))]);
        let original = source.clone();
        let mut edit = source.begin_batch_edit().expect("valid fixed topology");
        edit.remove_atom(AtomId::new(2)).expect("atom 2 exists");
        edit.remove_atom(AtomId::new(4)).expect("atom 4 exists");

        let (result, _) = edit.finish().expect("descending endpoint remaps succeed");

        assert_eq!(source, original);
        assert_eq!(
            result.bonds[0].prop("_MolFileBondEndPts"),
            Some(&crate::PropertyValue::String("(2 3 4)".to_owned()))
        );
        assert_eq!(
            result.bonds[0].prop("_MolFileBondAttach"),
            Some(&crate::PropertyValue::String("ANY".to_owned()))
        );
    }

    #[test]
    fn cf3d_frag_f13_zero_count_matching_endpoint_preserves_unsigned_decrement() {
        let source = f12_topology(5, vec![f13_bond(0, 0, 1, Some("(0 3 4)"), Some("ANY"))]);
        let original = source.clone();
        let mut edit = source.begin_batch_edit().expect("valid fixed topology");
        edit.remove_atom(AtomId::new(2)).expect("atom 2 exists");

        let (result, _) = edit.finish().expect("source unsigned decrement wraps");

        assert_eq!(source, original);
        assert_eq!(
            result.bonds[0].prop("_MolFileBondEndPts"),
            Some(&crate::PropertyValue::String("(4294967295 3)".to_owned()))
        );
        assert_eq!(
            result.bonds[0].prop("_MolFileBondAttach"),
            Some(&crate::PropertyValue::String("ANY".to_owned()))
        );
    }

    #[test]
    fn cf3d_frag_f13_errors_keep_source_bond_and_token_context_atomically() {
        let cases = [
            ("()", 0, "", BondEndPointsParseErrorKind::InvalidArgument),
            (
                "(abc 1e309)",
                0,
                "abc",
                BondEndPointsParseErrorKind::InvalidArgument,
            ),
            (
                "(1 abc)",
                1,
                "abc",
                BondEndPointsParseErrorKind::InvalidArgument,
            ),
            (
                "(1e309 1)",
                0,
                "1e309",
                BondEndPointsParseErrorKind::OutOfRange,
            ),
            (
                "(1 1e309)",
                1,
                "1e309",
                BondEndPointsParseErrorKind::OutOfRange,
            ),
            (
                "(-1 1)",
                0,
                "-1",
                BondEndPointsParseErrorKind::UnrepresentableEndpoint,
            ),
            (
                "(1 4294967296)",
                1,
                "4294967296",
                BondEndPointsParseErrorKind::UnrepresentableEndpoint,
            ),
        ];

        for (endpts, token_index, token, kind) in cases {
            let source = f12_topology(
                4,
                vec![
                    f13_bond(0, 0, 1, None, None),
                    f13_bond(1, 2, 3, Some(endpts), Some("ANY")),
                ],
            );
            let original = source.clone();
            let mut edit = source.begin_batch_edit().expect("valid fixed topology");
            edit.remove_atom(AtomId::new(0)).expect("atom 0 exists");

            let error = edit.finish().expect_err("bad retained ENDPTS fails");

            assert_eq!(source, original, "failure preserves source for {endpts:?}");
            assert_eq!(
                error,
                TopologyEditError::BondEndPointsParse {
                    bond: BondId::new(1),
                    token_index,
                    token: token.to_owned(),
                    kind,
                },
                "typed first-error context for {endpts:?}"
            );
        }
    }

    #[test]
    fn cf3d_frag_f14_batch_removal_remaps_sgroups_and_stereo_in_source_order() {
        let source = f14_topology();
        let original = source.clone();
        let mut edit = source.begin_batch_edit().expect("valid F14 topology");
        edit.remove_atom(AtomId::new(2))
            .expect("atom in the directly referenced group exists");
        edit.remove_bond(BondId::new(4))
            .expect("bond in the directly referenced group exists");

        let (result, mapping) = edit.finish().expect("group references remap");

        assert_eq!(source, original, "failed or successful edit is detached");
        assert_eq!(result.substance_groups.len(), 2);
        assert_eq!(
            result
                .substance_groups
                .iter()
                .map(SubstanceGroup::id)
                .collect::<Vec<_>>(),
            vec![SubstanceGroupId::new(0), SubstanceGroupId::new(1)]
        );

        let root = &result.substance_groups[0];
        assert_eq!(root.external_id(), Some(501));
        assert_eq!(root.rdkit_sequence_id(), Some(71));
        assert_eq!(root.label(), Some("retained-root"));
        assert_eq!(root.data_fields(), &["root-field"]);
        assert_eq!(root.props().get("vendor").map(String::as_str), Some("keep"));
        assert_eq!(
            root.atoms(),
            &[AtomId::new(6), AtomId::new(0), AtomId::new(4)]
        );
        assert_eq!(
            root.bonds(),
            &[BondId::new(2), BondId::new(0), BondId::new(1)]
        );
        assert_eq!(root.bond_role(BondId::new(2)), SGroupBondRole::Contained);
        assert_eq!(
            root.head_crossing_bonds(),
            &[BondId::new(2), BondId::new(1)]
        );
        assert_eq!(
            root.crossing_bond_correspondence(),
            &[BondId::new(1), BondId::new(0)]
        );
        assert_eq!(root.parent_atoms(), &[AtomId::new(5)]);
        assert_eq!(
            root.attach_points(),
            &[SGroupAttachPoint {
                atom: AtomId::new(4),
                leaving_atom: Some(AtomId::new(6)),
                label: Some("AP".into()),
                order: Some(4),
            }]
        );
        assert_eq!(
            root.cstates(),
            &[SGroupCState::new(BondId::new(1), [1.0, 2.0, 3.0])]
        );

        let retained_child = &result.substance_groups[1];
        assert_eq!(retained_child.parent(), Some(SubstanceGroupId::new(0)));
        assert_eq!(retained_child.atoms(), &[AtomId::new(5), AtomId::new(0)]);
        assert_eq!(retained_child.bonds(), &[BondId::new(2), BondId::new(1)]);

        assert_eq!(result.stereo_groups.len(), 2);
        assert_eq!(result.stereo_groups[0].id(), Some(41));
        assert_eq!(result.stereo_groups[0].write_id(), 9);
        assert_eq!(result.stereo_groups[0].kind(), StereoGroupKind::And);
        assert_eq!(
            result.stereo_groups[0].atoms(),
            &[AtomId::new(6), AtomId::new(0), AtomId::new(4)]
        );
        assert_eq!(
            result.stereo_groups[0].bonds(),
            &[BondId::new(2), BondId::new(0), BondId::new(1)]
        );
        assert_eq!(result.stereo_groups[1].id(), Some(42));
        assert_eq!(result.stereo_groups[1].write_id(), 12);
        assert_eq!(result.stereo_groups[1].kind(), StereoGroupKind::Or);
        assert_eq!(
            result.stereo_groups[1].atoms(),
            &[AtomId::new(4), AtomId::new(0)]
        );
        assert_eq!(
            result.stereo_groups[1].bonds(),
            &[BondId::new(1), BondId::new(0)]
        );
        assert_eq!(
            mapping.bonds.old_to_new,
            vec![
                Some(BondId::new(0)),
                None,
                Some(BondId::new(1)),
                Some(BondId::new(2)),
                None,
            ]
        );
        assert!(result.validate().is_ok());
    }

    #[test]
    fn cf3d_frag_f14_noop_preserves_group_order_and_both_stereo_ids() {
        let source = f14_topology();
        let (result, mapping) = source
            .begin_batch_edit()
            .expect("valid F14 topology")
            .finish()
            .expect("no-op batch edit preserves groups");

        assert_eq!(result, source);
        assert_eq!(
            mapping.atoms.old_to_new,
            (0..8)
                .map(|index| Some(AtomId::new(index)))
                .collect::<Vec<_>>()
        );
        assert_eq!(
            mapping.bonds.old_to_new,
            (0..5)
                .map(|index| Some(BondId::new(index)))
                .collect::<Vec<_>>()
        );
        assert_eq!(result.substance_groups.len(), 9);
        assert_eq!(result.stereo_groups.len(), 5);
        assert_eq!(result.stereo_groups[0].id(), Some(41));
        assert_eq!(result.stereo_groups[0].write_id(), 9);
        assert_eq!(result.stereo_groups[4].id(), Some(45));
        assert_eq!(result.stereo_groups[4].write_id(), 15);
    }
}
