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
    #[error("{0}")]
    StereoGroup(#[from] crate::StereoGroupError),

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
    #[error("bond {bond} stereo reference {atom} is out of range for {atom_count} atoms")]
    StereoReferenceOutOfRange {
        bond: BondId,
        atom: AtomId,
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
    #[error("{0}")]
    StereoGroup(#[from] crate::StereoGroupError),

    #[error("source coordinate edit failed: {0}")]
    Coordinate(#[from] crate::CoordinateValidationError),
    #[error("molecule property operation failed: {0}")]
    MoleculeProperty(#[from] crate::MoleculePropertyError),
    #[error("atom property operation failed: {0}")]
    AtomProperty(#[from] crate::AtomPropertyError),
    #[error("source SubstanceGroup {group:?} has a cyclic raw PARENT hierarchy")]
    SourceGroupParentCycle { group: crate::SubstanceGroupId },
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
        token: crate::PropertyText,
        kind: BondEndPointsParseErrorKind,
    },
    #[error("source batch atom and bond masks must both be active")]
    IncompleteSourceBatchMasks,
    #[error("NULL atom passed in")]
    NullBondAtom,
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
/// Original row counts establish the mapping. The owned working value and
/// pending masks remain detached; this type has no live Molecule authority.
#[derive(Debug, Clone)]
pub struct TopologyBatchEdit {
    source_atom_count: usize,
    source_bond_count: usize,
    working: TopologyBlock,
    added_neighbors: Vec<Vec<crate::NeighborRef>>,
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
        // The borrowed convenience entry still clones one working block to
        // isolate the caller. Original row counts require no second copy.
        self.clone().into_batch_edit()
    }

    /// Start the same detached editor by moving its sole working block.
    #[doc(hidden)]
    pub fn into_batch_edit(self) -> Result<TopologyBatchEdit, TopologyEditError> {
        // RDKit❗✔️: void RWMol::beginBatchEdit() {
        // RDKit❗✔️:   if (dp_delAtoms || dp_delBonds) {
        // RDKit❗✔️:     throw ValueErrorException("Attempt to re-enter batchEdit mode");
        // RDKit❗✔️:   }
        // RDKit❗✔️:   dp_delAtoms.reset(new boost::dynamic_bitset<>(getNumAtoms()));
        // RDKit❗✔️:   dp_delBonds.reset(new boost::dynamic_bitset<>(getNumBonds()));
        // RDKit❗✔️: }
        // An owned TopologyBlock is not already an editor, so re-entry cannot
        // occur. Linear validity/mask setup and moved ownership avoid cloning.
        self.validate().map_err(|error| match error {
            TopologyValidationError::StereoGroup(cause) => TopologyEditError::StereoGroup(cause),
            other => TopologyEditError::InvalidSource(other),
        })?;
        let source_atom_count = self.atoms.len();
        let source_bond_count = self.bonds.len();
        Ok(TopologyBatchEdit {
            source_atom_count,
            source_bond_count,
            remove_atoms: vec![false; source_atom_count],
            remove_bonds: vec![false; source_bond_count],
            added_neighbors: vec![Vec::new(); source_atom_count],
            working: self,
        })
    }

    pub fn reordered_atoms(
        &self,
        old_atom_order: &[AtomId],
    ) -> Result<(Self, TopologyMapping), TopologyEditError> {
        self.validate().map_err(|error| match error {
            TopologyValidationError::StereoGroup(cause) => TopologyEditError::StereoGroup(cause),
            other => TopologyEditError::InvalidSource(other),
        })?;
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
            .map(|group| group.remapped(&atom_old_to_new, &bond_old_to_new))
            .collect::<Result<Vec<_>, _>>()?
            .into_iter()
            .flatten()
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
            .map_err(|error| match error {
                TopologyValidationError::StereoGroup(cause) => {
                    TopologyEditError::StereoGroup(cause)
                }
                other => TopologyEditError::InvalidResult(other),
            })?;
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
            for atom in bond.stereo_atom_references() {
                if atom.index() >= self.atoms.len() {
                    return Err(TopologyValidationError::StereoReferenceOutOfRange {
                        bond: bond.id(),
                        atom: *atom,
                        atom_count: self.atoms.len(),
                    });
                }
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
            group.validate_members()?;
        }
        // Earlier checks preserve the original error precedence and establish
        // dense IDs, valid endpoints and no self-loops. Only the final CSR
        // comparison (including duplicate edges) is replaced, not skipped.
        if !self
            .adjacency
            .matches_validated_topology(self.atoms.len(), &self.bonds)
        {
            return Err(TopologyValidationError::AdjacencyMismatch);
        }
        Ok(())
    }
}

fn parse_endpts_decimal_prefix(token: &[u8]) -> Result<(f64, usize), BondEndPointsParseErrorKind> {
    // BEGIN RDKIT CPP FUNCTION RWMol::batchRemoveAtoms (ENDPTS conversions)
    // RDKit❗✔️: unsigned int num_ats = std::stod(*tokens.begin());
    // RDKit❗✔️: std::transform(beg, tokens.end(), std::back_inserter(oats),
    // RDKit❗✔️:                [](const std::string &a) { return std::stod(a); });
    // END RDKIT CPP FUNCTION
    // F13-NORMAL authorizes an atof-like decimal-prefix scope: ordinary finite
    // decimal values and reasonable scientific notation use Rust's binary64
    // conversion after this ASCII scan. This does not claim equivalence for
    // every std::stod/libc spelling or rare rounding boundary.
    // A scanned normal token is O(n). Only its grammar-selected ASCII number
    // is copied to the same Rust float primitive; arbitrary property bytes
    // are never decoded or rejected as text. The existing F13 numeric scope
    // and framing remain unchanged; this adds an O(prefix) primitive buffer.
    let bytes = token;
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
    let numeric: String = token[number_start..index]
        .iter()
        .map(|&byte| char::from(byte))
        .collect();
    let parsed = numeric
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
    /// Move the actual detached graph and active source deletion masks.
    /// This value projection performs no commit, clone, or metadata invention.
    #[doc(hidden)]
    pub fn into_source_batch_parts(self) -> (TopologyBlock, Option<Vec<bool>>, Option<Vec<bool>>) {
        (
            self.working,
            Some(self.remove_atoms),
            Some(self.remove_bonds),
        )
    }

    /// Borrow one checked atom row of this detached working value.
    #[doc(hidden)]
    pub fn atom_mut(&mut self, atom: AtomId) -> Result<&mut Atom, TopologyEditError> {
        // RDKit❗✔️: Atom *ROMol::getAtomWithIdx(unsigned int idx) {
        // RDKit❗✔️:   URANGE_CHECK(idx, getNumAtoms());
        // RDKit❗✔️:
        // RDKit❗✔️:   auto vd = boost::vertex(idx, d_graph);
        // RDKit❗✔️:   auto res = d_graph[vd];
        // RDKit❗✔️:   POSTCONDITION(res, "");
        // RDKit❗✔️:   return res;
        // RDKit❗✔️: }
        // RDKit❗✔️:
        // RDKit❗✔️: const Atom *ROMol::getAtomWithIdx(unsigned int idx) const {
        // RDKit❗✔️:   URANGE_CHECK(idx, getNumAtoms());
        let atom_count = self.working.atoms.len();
        self.working
            .atoms
            .get_mut(atom.index())
            .ok_or(TopologyEditError::AtomOutOfRange { atom, atom_count })
    }

    /// Borrow one checked bond row without exposing the editor's storage.
    #[doc(hidden)]
    pub fn bond_mut(&mut self, bond: BondId) -> Result<&mut Bond, TopologyEditError> {
        // RDKit❗✔️: Bond *ROMol::getBondWithIdx(unsigned int idx) {
        // RDKit❗✔️:   return const_cast<Bond *>(static_cast<const ROMol *>(this)->getBondWithIdx(
        // RDKit❗✔️:       idx));  // avoid code duplication
        // RDKit❗✔️: }
        // RDKit❗✔️:
        // RDKit❗✔️: const Bond *ROMol::getBondWithIdx(unsigned int idx) const {
        // RDKit❗✔️:   URANGE_CHECK(idx, getNumBonds());
        // RDKit❗✔️:
        // RDKit❗✔️:   // boost::graph doesn't give us random-access to edges,
        // RDKit❗✔️:   // so we have to iterate to it
        // RDKit❗✔️:   auto [iter, end] = getEdges();
        // RDKit❗✔️:   for (unsigned int i = 0; i < idx; i++) {
        // RDKit❗✔️:     ++iter;
        // RDKit❗✔️:   }
        // RDKit❗✔️:   const Bond *res = d_graph[*iter];
        // RDKit❗✔️:
        // RDKit❗✔️:   POSTCONDITION(res != nullptr, "Invalid bond requested");
        // RDKit❗✔️:   return res;
        // RDKit❗✔️: }
        // Indexed canonical rows avoid ROMol's linear edge-index walk.
        let bond_count = self.working.bonds.len();
        self.working
            .bonds
            .get_mut(bond.index())
            .ok_or(TopologyEditError::BondOutOfRange { bond, bond_count })
    }

    /// Replace one detached bond using the canonical algorithm readers supplied
    /// by its consumer. Same stable ID also preserves bookmark references.
    #[allow(dead_code)]
    pub(crate) fn replace_bond_source<E: From<TopologyEditError>>(
        &mut self,
        bond: BondId,
        replacement: &Bond,
        preserve_props: bool,
        keep_sgroups: bool,
        order_reader: impl FnMut(crate::BondOrder) -> Result<f64, E>,
        uint_reader: impl FnMut(&crate::PropertyValue) -> Result<u32, E>,
    ) -> Result<(), E> {
        replace_source_bond(
            &mut self.working,
            bond,
            replacement,
            preserve_props,
            keep_sgroups,
            order_reader,
            uint_reader,
        )
    }

    /// Source-visible edges include additions and pending removals until finish.
    #[doc(hidden)]
    pub fn bond_between_atoms(
        &self,
        begin: AtomId,
        end: AtomId,
    ) -> Result<Option<BondId>, TopologyEditError> {
        source_bond_between_atoms(self.working.atoms.len(), begin, end, || {
            self.working
                .adjacency
                .neighbors_of(begin.index())
                .iter()
                .chain(self.added_neighbors[begin.index()].iter())
        })
    }

    // Private source-overload dependency with explicit detached side state.
    // No new facade or operation capability is exposed by this implementation.
    fn add_default_atom_with_source_state(
        &mut self,
        update_label: bool,
        bookmarks: &mut std::collections::BTreeMap<i32, Vec<AtomId>>,
        coordinates: &mut crate::CoordinateBlock,
    ) -> Result<AtomId, crate::CoordinateValidationError> {
        // BEGIN RDKIT CPP FUNCTION RWMol::addAtom(bool)
        // RDKit❗✔️: unsigned int RWMol::addAtom(bool updateLabel) {
        // RDKit❗✔️:   auto *atom_p = new Atom();
        // RDKit❗✔️:   atom_p->setOwningMol(this);
        // RDKit❗✔️:   auto which = boost::add_vertex(d_graph);
        // RDKit❗✔️:   d_graph[which] = atom_p;
        // RDKit❗✔️:   atom_p->setIdx(which);
        // RDKit❗✔️:   if (updateLabel) {
        // RDKit❗✔️:     clearAtomBookmark(ci_RIGHTMOST_ATOM);
        // RDKit❗✔️:     setAtomBookmark(atom_p, ci_RIGHTMOST_ATOM);
        // RDKit❗✔️:   }
        // RDKit❗✔️:
        // RDKit❗✔️:   // add atom to any conformers as well, if we have any
        // RDKit❗✔️:   for (auto &conf : d_confs) {
        // RDKit❗✔️:     conf->setAtomPos(which, RDGeom::Point3D(0.0, 0.0, 0.0));
        // RDKit❗✔️:   }
        // RDKit❗✔️:   return rdcast<unsigned int>(which);
        // RDKit❗✔️: }
        // END RDKIT CPP FUNCTION RWMol::addAtom(bool)
        // BEGIN RDKIT REACHED Atom default constructor/initAtom
        // RDKit❗✔️: Atom::Atom() : RDProps() {
        // RDKit❗✔️:   d_atomicNum = 0;
        // RDKit❗✔️:   initAtom();
        // RDKit❗✔️: }
        // RDKit❗✔️: void Atom::initAtom() {
        // RDKit❗✔️:   df_isAromatic = false;
        // RDKit❗✔️:   df_noImplicit = false;
        // RDKit❗✔️:   d_numExplicitHs = 0;
        // RDKit❗✔️:   d_numRadicalElectrons = 0;
        // RDKit❗✔️:   d_formalCharge = 0;
        // RDKit❗✔️:   d_index = 0;
        // RDKit❗✔️:   d_isotope = 0;
        // RDKit❗✔️:   d_chiralTag = CHI_UNSPECIFIED;
        // RDKit❗✔️:   d_hybrid = UNSPECIFIED;
        // RDKit❗✔️:   dp_mol = nullptr;
        // RDKit❗✔️:   dp_monomerInfo = nullptr;
        // RDKit❗✔️:
        // RDKit❗✔️:   d_implicitValence = -1;
        // RDKit❗✔️:   d_explicitValence = -1;
        // RDKit❗✔️: }
        // RDKit❗✔️:
        // RDKit❗✔️: Atom::~Atom() { delete dp_monomerInfo; }
        // END RDKIT REACHED Atom default constructor/initAtom
        // Atom::Atom's zero atomic number is the existing Element::DUMMY,
        // not a guessed element. AtomSpec's canonical constructor supplies the
        // native initAtom fields; dense append uses the sole topology owner.
        // Bookmark clear then append happens after graph insertion, before any
        // conformer write. Later coordinate errors retain earlier detached
        // effects, like native mutation; runtime commit authority stays outside.
        // Detached stable IDs replace owning pointers; native uint32 graph
        // width/debug casts remain a boundary gap. Dimension rows commute for
        // these independent zero writes (overflow fails before any vector write),
        // so splitting their storage does not alter modeled successful state.
        // Cost: amortized graph append, ordered-map bookmark erase/insert and
        // amortized append to each coordinate Vec, like native graph/map/list.
        let atom = self.add_atom(AtomSpec::new(cosmolkit_types::Element::DUMMY));
        if update_label {
            const RIGHTMOST_ATOM: i32 = -0xBADBEEF;
            // RDKit❗✔️: void ROMol::clearAtomBookmark(const int mark) { d_atomBookmarks.erase(mark); }
            bookmarks.remove(&RIGHTMOST_ATOM);
            // RDKit❗✔️:   void setAtomBookmark(Atom *at, int mark) {
            // RDKit❗✔️:     d_atomBookmarks[mark].push_back(at);
            // RDKit❗✔️:   }
            bookmarks.entry(RIGHTMOST_ATOM).or_default().push(atom);
        }
        for conformer in &mut coordinates.conformers_2d {
            conformer.source_set_atom_position(atom.index(), [0.0; 2])?;
        }
        for conformer in &mut coordinates.conformers_3d {
            conformer.source_set_atom_position(atom.index(), [0.0; 3])?;
        }
        Ok(atom)
    }

    pub fn add_atom(&mut self, spec: AtomSpec) -> AtomId {
        // RDKit✔️✔️: if (dp_delAtoms->size() < getNumAtoms()) {
        // RDKit✔️✔️:   dp_delAtoms->resize(getNumAtoms());
        // RDKit✔️✔️: }
        let id = AtomId::new(self.working.atoms.len());
        self.working.atoms.push(Atom::from_spec(id, spec));
        self.remove_atoms.push(false);
        self.added_neighbors.push(Vec::new());
        id
    }

    #[doc(hidden)]
    pub fn add_bond_from_atoms(
        &mut self,
        begin: Option<&Atom>,
        end: Option<&Atom>,
        order: crate::BondOrder,
    ) -> Result<usize, TopologyEditError> {
        // RDKit❗✔️: unsigned int RWMol::addBond(Atom *atom1, Atom *atom2, Bond::BondType bondType) {
        // RDKit❗✔️:   PRECONDITION(atom1 && atom2, "NULL atom passed in");
        // RDKit❗✔️:   return addBond(atom1->getIdx(), atom2->getIdx(), bondType);
        // RDKit❗✔️: }
        // Null errors precede index lookup. Foreign atoms contribute only their
        // IDs, as the source has no ownership predicate in this overload.
        // Cost is constant two field reads and sole-kernel delegation; all
        // native pointer/index and order-overload boundaries stay explicit.
        let (Some(begin), Some(end)) = (begin, end) else {
            return Err(TopologyEditError::NullBondAtom);
        };
        self.add_bond_by_order(begin.id(), end.id(), order)
    }

    #[doc(hidden)]
    pub fn add_bond_by_order(
        &mut self,
        begin: AtomId,
        end: AtomId,
        order: crate::BondOrder,
    ) -> Result<usize, TopologyEditError> {
        let count = add_source_bond_order(
            &mut self.working.atoms,
            &mut self.working.bonds,
            SourceBondNeighbors {
                original: Some(&self.working.adjacency),
                appended: &mut self.added_neighbors,
            },
            Some(SourceBondBatchMasks {
                atoms: &self.remove_atoms,
                bonds: &mut self.remove_bonds,
            }),
            begin,
            end,
            order,
        )?;
        // Editor masks are always aligned, unlike native lazily resized masks.
        // This canonical editor invariant adds only zero slots after the source
        // helper and cannot change an already scheduled deletion.
        self.remove_bonds.resize(count, false);
        Ok(count)
    }

    /// Source pointer overload for a detached bond with source constructor defaults.
    /// Delegate the sole value kernel without rich-spec or aromatic projection.
    #[doc(hidden)]
    pub fn add_bond_value_source(
        &mut self,
        bond: std::borrow::Cow<'_, Bond>,
    ) -> Result<usize, TopologyEditError> {
        let count = add_source_bond_value(
            self.working.atoms.len(),
            &mut self.working.bonds,
            SourceBondNeighbors {
                original: Some(&self.working.adjacency),
                appended: &mut self.added_neighbors,
            },
            bond,
        )?;
        // Editor deletion masks are aligned side state. Native pointer add does
        // not inherit endpoint atom deletions, unlike the order overload.
        self.remove_bonds.resize(count, false);
        Ok(count)
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
        // This editor additionally validates rich detached specs above and
        // projects ordinary aromatic flags below. The canonical pointer/value
        // overload itself does not alter atoms, valence caches or batch masks.
        let begin = spec.begin().index();
        let end = spec.end().index();
        let aromatic_order = spec.query().is_none() && spec.order() == crate::BondOrder::Aromatic;
        let count = add_source_bond_value(
            self.working.atoms.len(),
            &mut self.working.bonds,
            SourceBondNeighbors {
                original: Some(&self.working.adjacency),
                appended: &mut self.added_neighbors,
            },
            std::borrow::Cow::Owned(Bond::from_spec(BondId::new(0), spec)),
        )?;
        let id = BondId::new(count - 1);
        // These infallible reflected fields have no intervening error boundary.
        if aromatic_order {
            self.working.bonds[id.index()].set_aromatic(true);
            self.working.atoms[begin].set_aromatic(true);
            self.working.atoms[end].set_aromatic(true);
        }
        self.remove_bonds
            .push(self.remove_atoms[begin] || self.remove_atoms[end]);
        Ok(id)
    }

    // Project source graph additions into the detached CSR before delegating
    // an algorithm that reads topology adjacency directly. No logical edge or
    // pending deletion changes here; avoid rebuilding when no additions exist.
    fn synchronize_source_neighbors(&mut self) {
        if self.added_neighbors.iter().any(|row| !row.is_empty()) {
            self.working.adjacency =
                AdjacencyList::from_topology(self.working.atoms.len(), &self.working.bonds);
            for row in &mut self.added_neighbors {
                row.clear();
            }
        }
    }

    pub fn remove_atom(&mut self, atom: AtomId) -> Result<(), TopologyEditError> {
        if atom.index() < self.working.atoms.len() {
            self.synchronize_source_neighbors();
        }
        remove_atom_source_index(
            &mut self.working,
            atom,
            SourceAtomRemovalState::Batch {
                atoms: &mut self.remove_atoms,
                bonds: &mut self.remove_bonds,
            },
        )
    }

    #[doc(hidden)]
    pub fn remove_bond_between_atoms(
        &mut self,
        begin: AtomId,
        end: AtomId,
    ) -> Result<(), TopologyEditError> {
        if begin.index() < self.working.atoms.len() && end.index() < self.working.atoms.len() {
            self.synchronize_source_neighbors();
        }
        remove_bond_source(
            &mut self.working,
            begin,
            end,
            SourceBondRemovalState::Batch {
                pending: &mut self.remove_bonds,
            },
        )
    }

    pub fn remove_bond(&mut self, bond: BondId) -> Result<(), TopologyEditError> {
        let bond_count = self.working.bonds.len();
        let Some(slot) = self.remove_bonds.get_mut(bond.index()) else {
            return Err(TopologyEditError::BondOutOfRange { bond, bond_count });
        };
        *slot = true;
        Ok(())
    }

    /// Finish using the sole source batch-commit kernel. Explicit detached
    /// side values carry source coordinate/property errors without live authority.
    #[doc(hidden)]
    pub fn finish_source<E: From<TopologyEditError>>(
        mut self,
        coordinates: &mut crate::CoordinateBlock,
        properties: &mut crate::MoleculeProperties,
        uint_reader: &mut dyn FnMut(&crate::PropertyValue) -> Result<u32, E>,
        text_reader: &mut dyn FnMut(&crate::PropertyValue) -> Result<crate::PropertyText, E>,
    ) -> Result<(TopologyBlock, TopologyMapping), E> {
        self.synchronize_source_neighbors();
        // Structural row transport only: source deletion leaves retained rows
        // in encounter order, with appended rows lacking an old counterpart.
        let mut mapping = TopologyMapping {
            atoms: AtomMapping {
                old_to_new: vec![None; self.source_atom_count],
                new_to_old: Vec::new(),
            },
            bonds: BondMapping {
                old_to_new: vec![None; self.source_bond_count],
                new_to_old: Vec::new(),
            },
        };
        for (old, removed) in self.remove_atoms.iter().copied().enumerate() {
            if removed {
                continue;
            }
            let new = AtomId::new(mapping.atoms.new_to_old.len());
            let source = (old < self.source_atom_count).then(|| AtomId::new(old));
            if let Some(old) = source {
                mapping.atoms.old_to_new[old.index()] = Some(new);
            }
            mapping.atoms.new_to_old.push(source);
        }
        for (old, removed) in self.remove_bonds.iter().copied().enumerate() {
            if removed {
                continue;
            }
            let new = BondId::new(mapping.bonds.new_to_old.len());
            let source = (old < self.source_bond_count).then(|| BondId::new(old));
            if let Some(old) = source {
                mapping.bonds.old_to_new[old.index()] = Some(new);
            }
            mapping.bonds.new_to_old.push(source);
        }
        let mut atoms = Some(self.remove_atoms);
        let mut bonds = Some(self.remove_bonds);
        commit_batch_edit_source(
            &mut self.working,
            SourceBatchCommitState {
                atoms: &mut atoms,
                bonds: &mut bonds,
                atom_bookmarks: None,
                bond_bookmarks: None,
                coordinates,
                properties,
                reset_ring: None,
                uint_reader,
                text_reader,
            },
        )?;
        self.working
            .validate()
            .map_err(|e| E::from(TopologyEditError::InvalidResult(e)))?;
        mapping
            .validate_for_counts(
                self.source_atom_count,
                self.working.atoms.len(),
                self.source_bond_count,
                self.working.bonds.len(),
            )
            .map_err(|e| E::from(TopologyEditError::InvalidMapping(e)))?;
        Ok((self.working, mapping))
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
        // BEGIN RDKIT CPP FUNCTION: third_party/rdkit/Code/RDGeneral/RDValue.h:192-266
        // RDKit❗✔️: inline bool rdvalue_tostring(RDValue_cast_t val, std::string &res) {
        // RDKit❗✔️:   switch (val.getTag()) {
        // RDKit❗✔️:     case RDTypeTag::StringTag:
        // RDKit❗✔️:       res = rdvalue_cast<std::string>(val);
        // RDKit❗✔️:       break;
        // RDKit❗✔️:     case RDTypeTag::IntTag:
        // RDKit❗✔️:       res = boost::lexical_cast<std::string>(rdvalue_cast<int>(val));
        // RDKit❗✔️:       break;
        // RDKit❗✔️:     case RDTypeTag::DoubleTag: {
        // RDKit❗✔️:       Utils::LocaleSwitcher ls;  // for lexical cast...
        // RDKit❗✔️:       res = boost::lexical_cast<std::string>(rdvalue_cast<double>(val));
        // RDKit❗✔️:       break;
        // RDKit❗✔️:     }
        // RDKit❗✔️:     case RDTypeTag::UnsignedIntTag:
        // RDKit❗✔️:       res = boost::lexical_cast<std::string>(rdvalue_cast<unsigned int>(val));
        // RDKit❗✔️:       break;
        // RDKit❗✔️: #ifdef RDVALUE_HASBOOL
        // RDKit❗✔️:     case RDTypeTag::BoolTag:
        // RDKit❗✔️:       res = boost::lexical_cast<std::string>(rdvalue_cast<bool>(val));
        // RDKit❗✔️:       break;
        // RDKit❗✔️: #endif
        // RDKit❗✔️:     case RDTypeTag::FloatTag: {
        // RDKit❗✔️:       Utils::LocaleSwitcher ls;  // for lexical cast...
        // RDKit❗✔️:       res = boost::lexical_cast<std::string>(rdvalue_cast<float>(val));
        // RDKit❗✔️:       break;
        // RDKit❗✔️:     }
        // RDKit❗✔️:     case RDTypeTag::VecDoubleTag: {
        // RDKit❗✔️:       // vectToString uses std::imbue for locale
        // RDKit❗✔️:       res = vectToString<double>(val);
        // RDKit❗✔️:       break;
        // RDKit❗✔️:     }
        // RDKit❗✔️:     case RDTypeTag::VecFloatTag: {
        // RDKit❗✔️:       // vectToString uses std::imbue for locale
        // RDKit❗✔️:       res = vectToString<float>(val);
        // RDKit❗✔️:       break;
        // RDKit❗✔️:     }
        // RDKit❗✔️:     case RDTypeTag::VecIntTag:
        // RDKit❗✔️:       res = vectToString<int>(val);
        // RDKit❗✔️:       break;
        // RDKit❗✔️:     case RDTypeTag::VecUnsignedIntTag:
        // RDKit❗✔️:       res = vectToString<unsigned int>(val);
        // RDKit❗✔️:       break;
        // RDKit❗✔️:     case RDTypeTag::VecStringTag:
        // RDKit❗✔️:       res = vectToString<std::string>(val);
        // RDKit❗✔️:       break;
        // RDKit❗✔️:     case RDTypeTag::AnyTag: {
        // RDKit❗✔️:       Utils::LocaleSwitcher ls;  // for lexical cast...
        // RDKit❗✔️:       try {
        // RDKit❗✔️:         res = std::any_cast<std::string>(rdvalue_cast<std::any &>(val));
        // RDKit❗✔️:       } catch (const std::bad_any_cast &) {
        // RDKit❗✔️:         auto &rdtype = rdvalue_cast<std::any &>(val).type();
        // RDKit❗✔️:         if (rdtype == typeid(long)) {
        // RDKit❗✔️:           res = boost::lexical_cast<std::string>(
        // RDKit❗✔️:               std::any_cast<long>(rdvalue_cast<std::any &>(val)));
        // RDKit❗✔️:         } else if (rdtype == typeid(int64_t)) {
        // RDKit❗✔️:           res = boost::lexical_cast<std::string>(
        // RDKit❗✔️:               std::any_cast<int64_t>(rdvalue_cast<std::any &>(val)));
        // RDKit❗✔️:         } else if (rdtype == typeid(uint64_t)) {
        // RDKit❗✔️:           res = boost::lexical_cast<std::string>(
        // RDKit❗✔️:               std::any_cast<uint64_t>(rdvalue_cast<std::any &>(val)));
        // RDKit❗✔️:         } else if (rdtype == typeid(unsigned long)) {
        // RDKit❗✔️:           res = boost::lexical_cast<std::string>(
        // RDKit❗✔️:               std::any_cast<unsigned long>(rdvalue_cast<std::any &>(val)));
        // RDKit❗✔️:         } else {
        // RDKit❗✔️:           throw;
        // RDKit❗✔️:           return false;
        // RDKit❗✔️:         }
        // RDKit❗✔️:       }
        // RDKit❗✔️:       break;
        // RDKit❗✔️:     }
        // RDKit❗✔️:     default:
        // RDKit❗✔️:       res = "";
        // RDKit❗✔️:   }
        // RDKit❗✔️:   return true;
        // RDKit❗✔️: }
        // END RDKIT CPP FUNCTION: third_party/rdkit/Code/RDGeneral/RDValue.h:192-266
        // UInt formats as unsigned decimal without parentheses. Preserve the
        // source no-read ENDPTS guard without allocating its unused text.
        for removed_atom_index in (0..self.remove_atoms.len()).rev() {
            if !self.remove_atoms[removed_atom_index] {
                continue;
            }
            for (bond_index, bond) in self.working.bonds.iter_mut().enumerate() {
                if self.remove_bonds[bond_index] {
                    continue;
                }
                let Some(source_property) = bond.prop("_MolFileBondEndPts") else {
                    continue;
                };
                let source_property = match source_property {
                    crate::PropertyValue::String(value) => value.clone(),
                    // RDKit✔️✔️: rdvalue_tostring(i.val, res);
                    // RDKit✔️✔️: if ('(' == sprop.front() && ')' == sprop.back()) {
                    // Dict::getValIfPresent(string&) projects Int/UInt/Double/Bool
                    // through lexical_cast. None of those scalar spellings
                    // begins with '(' (including signed zero, Inf and NaN),
                    // so the source skips this branch. Avoid that unobserved
                    // allocation without moving formatting into the model.
                    // String-vector formatting begins with [, so the source
                    // parenthesis guard also skips it for every raw payload.
                    crate::PropertyValue::StringVector(_)
                    | crate::PropertyValue::IntVector(_)
                    | crate::PropertyValue::Int(_)
                    | crate::PropertyValue::UInt(_)
                    | crate::PropertyValue::Double(_)
                    | crate::PropertyValue::Bool(_) => continue,
                };
                update_source_bond_endpoints_property(bond, removed_atom_index, &source_property)?;
            }
        }
        let old_atom_count = self.source_atom_count;
        let old_bond_count = self.source_bond_count;
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
            .map(|source_group| {
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
                if group.is_empty() {
                    Ok(None)
                } else {
                    group.remapped(&all_atom_to_new, &bond_old_to_new)
                }
            })
            .collect::<Result<Vec<_>, _>>()?
            .into_iter()
            .flatten()
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
        let mut topology =
            TopologyBlock::try_from_parts(atoms, bonds, substance_groups, stereo_groups).map_err(
                |error| match error {
                    TopologyValidationError::StereoGroup(cause) => {
                        TopologyEditError::StereoGroup(cause)
                    }
                    other => TopologyEditError::InvalidResult(other),
                },
            )?;
        // RDKit❗✔️:   // fix properties
        // RDKit❗✔️:   clearComputedProps(true);
        // RDKit❗✔️:   for (auto atom : atoms()) {
        // RDKit❗✔️:     atom->clearComputedProps();
        // RDKit❗✔️:   }
        // RDKit❗✔️:   for (auto bond : bonds()) {
        // RDKit❗✔️:     bond->clearComputedProps();
        // RDKit❗✔️:   }
        // RWMol's no-removal branch preserves these properties; molecule
        // properties and ring validity stay at their canonical runtime owner.
        if self.remove_atoms.iter().any(|removed| *removed)
            || self.remove_bonds.iter().any(|removed| *removed)
        {
            for atom in &mut topology.atoms {
                atom.clear_computed_props()?;
            }
            for bond in &mut topology.bonds {
                bond.clear_computed_props()
                    .map_err(TopologyEditError::InvalidBond)?;
            }
        }
        Ok((topology, mapping))
    }
}

#[cfg(test)]
mod tests {
    fn fixture_text(value: &crate::PropertyText) -> &str {
        std::str::from_utf8(value.as_bytes()).expect("unchanged UTF-8 fixture bytes")
    }

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
        root.set_prop("vendor", "keep")
            .expect("original fixture SGroup property write succeeds");

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
            .expect("valid distinct stereo members")
            .with_id(41)
            .with_write_id(9),
            crate::StereoGroup::new(
                StereoGroupKind::Or,
                vec![AtomId::new(5), AtomId::new(2), AtomId::new(0)],
                vec![BondId::new(2), BondId::new(4), BondId::new(0)],
            )
            .expect("valid distinct stereo members")
            .with_id(42)
            .with_write_id(12),
            crate::StereoGroup::new(StereoGroupKind::Absolute, vec![AtomId::new(2)], vec![])
                .expect("valid distinct stereo members")
                .with_id(43)
                .with_write_id(13),
            crate::StereoGroup::new(StereoGroupKind::Or, vec![], vec![BondId::new(4)])
                .expect("valid distinct stereo members")
                .with_id(44)
                .with_write_id(14),
            crate::StereoGroup::new(
                StereoGroupKind::And,
                vec![AtomId::new(2)],
                vec![BondId::new(4)],
            )
            .expect("valid distinct stereo members")
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
    fn direct_adjacency_validation_preserves_topology_error_precedence() {
        let mut block = f12_topology(
            3,
            vec![
                f12_bond(0, 0, 1, BondOrder::Single, BondStereo::None, None),
                f12_bond(1, 1, 2, BondOrder::Single, BondStereo::None, None),
            ],
        );
        // Duplicate edges are still an adjacency error, including reversed
        // endpoints. Earlier ID/endpoint failures must not be masked by it.
        block.bonds[1] = f12_bond(1, 1, 0, BondOrder::Single, BondStereo::None, None);
        assert_eq!(
            block.validate(),
            Err(TopologyValidationError::AdjacencyMismatch)
        );
        block.bonds[1] = f12_bond(0, 1, 0, BondOrder::Single, BondStereo::None, None);
        assert_eq!(
            block.validate(),
            Err(TopologyValidationError::BondIdMismatch {
                position: 1,
                id: BondId::new(0)
            })
        );
        block.bonds[1] = f12_bond(1, 1, 3, BondOrder::Single, BondStereo::None, None);
        assert_eq!(
            block.validate(),
            Err(TopologyValidationError::BondEndpointOutOfRange {
                bond: BondId::new(1),
                endpoint: "end",
                atom: AtomId::new(3),
                atom_count: 3,
            })
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
            let (value, consumed) = parse_endpts_decimal_prefix(token.as_bytes())
                .expect("fixed decimal prefix converts");
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
            let (value, consumed) = parse_endpts_decimal_prefix(token.as_bytes())
                .expect("C whitespace precedes a number");
            assert_eq!(value.to_bits(), 0x3ff0_0000_0000_0000, "{token:?}");
            assert_eq!(consumed, expected_consumed, "{token:?}");
        }

        for token in ["", " ", "\t\n\x0b\x0c\r", "abc", "+", "-", ".", "e1", "\0"] {
            assert_eq!(
                parse_endpts_decimal_prefix(token.as_bytes()),
                Err(BondEndPointsParseErrorKind::InvalidArgument),
                "no-conversion classification for {token:?}"
            );
        }
        assert_eq!(
            parse_endpts_decimal_prefix("1e309".as_bytes()),
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
        let (hex_prefix, hex_consumed) = parse_endpts_decimal_prefix("0x1p+2".as_bytes()).unwrap();
        assert_eq!(hex_prefix.to_bits(), 0x0000_0000_0000_0000);
        assert_eq!(hex_consumed, 1);

        let (rust_subnormal, subnormal_consumed) =
            parse_endpts_decimal_prefix("1e-310".as_bytes()).unwrap();
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
            Some(&crate::PropertyValue::String("(2 3 4)".into()))
        );
        assert_eq!(
            result.bonds[0].prop("_MolFileBondAttach"),
            Some(&crate::PropertyValue::String("ANY".into()))
        );
        assert_eq!(result.bonds[1].prop("_MolFileBondEndPts"), None);
        assert_eq!(result.bonds[1].prop("_MolFileBondAttach"), None);
        assert_eq!(
            result.bonds[2].prop("_MolFileBondEndPts"),
            Some(&crate::PropertyValue::String("raw 4 5".into()))
        );
        assert_eq!(
            result.bonds[2].prop("_MolFileBondAttach"),
            Some(&crate::PropertyValue::String("ANY".into()))
        );
        assert_eq!(
            result.bonds[3].prop("_MolFileBondEndPts"),
            Some(&crate::PropertyValue::String("(2 2 3)".into()))
        );
        assert_eq!(
            result.bonds[3].prop("_MolFileBondAttach"),
            Some(&crate::PropertyValue::String("ANY".into()))
        );
        assert_eq!(
            result.bonds[4].prop("_MolFileBondEndPts"),
            Some(&crate::PropertyValue::String("(8 7)".into()))
        );
        assert_eq!(
            result.bonds[4].prop("_MolFileBondAttach"),
            Some(&crate::PropertyValue::String("ANY".into()))
        );
    }

    #[test]
    fn uff_sync_endpts_typed_scalars_skip_rewrite_and_keep_exact_values() {
        let scalar_values = [
            (
                "nonparenthesized String",
                crate::PropertyValue::String("raw 4 5".into()),
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
                Some(&crate::PropertyValue::String("ANY".into())),
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
            Some(&crate::PropertyValue::String("(2 1 2)".into()))
        );
        assert_eq!(
            result.bonds[0].prop("_MolFileBondAttach"),
            Some(&crate::PropertyValue::String("ANY".into()))
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
            Some(&crate::PropertyValue::String("(2 3 4)".into()))
        );
        assert_eq!(
            result.bonds[0].prop("_MolFileBondAttach"),
            Some(&crate::PropertyValue::String("ANY".into()))
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
            Some(&crate::PropertyValue::String("(4294967295 3)".into()))
        );
        assert_eq!(
            result.bonds[0].prop("_MolFileBondAttach"),
            Some(&crate::PropertyValue::String("ANY".into()))
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
                    token: token.into(),
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
        assert_eq!(root.label().map(fixture_text), Some("retained-root"));
        assert_eq!(
            root.data_fields()
                .iter()
                .map(fixture_text)
                .collect::<Vec<_>>(),
            &["root-field"]
        );
        assert_eq!(
            root.props()
                .get("vendor".as_bytes())
                .map(|value| fixture_text(value.as_string().unwrap())),
            Some("keep")
        );
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

#[cfg(test)]
mod source_default_add_atom_complete_tests {
    use super::*;
    use crate::{Conformer2D, Conformer3D, CoordinateBlock, CoordinateDimension};
    use cosmolkit_types::Element;
    use std::collections::BTreeMap;

    #[test]
    fn default_atom_updates_label_and_zero_extends_every_conformer_without_changing_metadata() {
        let source = TopologyBlock::try_from_parts(
            vec![Atom::from_spec(AtomId::new(0), AtomSpec::new(Element::C))],
            vec![],
            vec![],
            vec![],
        )
        .unwrap();
        let mut edit = source.begin_batch_edit().unwrap();
        let mut bookmarks = BTreeMap::from([
            (-0xBADBEEF, vec![AtomId::new(0), AtomId::new(0)]),
            (17, vec![AtomId::new(0)]),
        ]);
        let mut coordinates = CoordinateBlock {
            conformers_2d: vec![Conformer2D::new(7, vec![[-0.0, 3.0]]).with_prop("note", "2D")],
            conformers_3d: vec![
                Conformer3D::new(7, vec![[1.0, 2.0, 3.0]], false).with_prop("note", "3D"),
            ],
            source_coordinate_dim: Some(CoordinateDimension::ThreeD),
            source_conformer_order: Some(vec![
                CoordinateDimension::ThreeD,
                CoordinateDimension::TwoD,
            ]),
        };
        let original_order = coordinates.source_conformer_order.clone();
        let atom = edit
            .add_default_atom_with_source_state(true, &mut bookmarks, &mut coordinates)
            .unwrap();
        assert_eq!(atom, AtomId::new(1));
        assert_eq!(
            edit.working.atoms[1],
            Atom::from_spec(atom, AtomSpec::new(Element::DUMMY))
        );
        assert_eq!(edit.working.atoms[0], source.atoms[0]);
        assert_eq!(bookmarks.get(&-0xBADBEEF), Some(&vec![atom]));
        assert_eq!(bookmarks.get(&17), Some(&vec![AtomId::new(0)]));
        assert_eq!(
            coordinates.conformers_2d[0].coordinates(),
            [[-0.0, 3.0], [0.0, 0.0]]
        );
        assert_eq!(
            coordinates.conformers_2d[0].coordinates()[0][0].to_bits(),
            (-0.0_f64).to_bits()
        );
        assert_eq!(
            coordinates.conformers_3d[0].coordinates(),
            [[1.0, 2.0, 3.0], [0.0, 0.0, 0.0]]
        );
        assert_eq!(coordinates.conformers_2d[0].id(), 7);
        assert_eq!(coordinates.conformers_3d[0].id(), 7);
        assert!(!coordinates.conformers_3d[0].is_3d());
        assert_eq!(
            coordinates.conformers_2d[0].props().get(b"note".as_slice()),
            Some(&crate::PropertyText::from("2D"))
        );
        assert_eq!(
            coordinates.conformers_3d[0].props().get(b"note".as_slice()),
            Some(&crate::PropertyText::from("3D"))
        );
        assert_eq!(coordinates.source_conformer_order, original_order);
        assert_eq!(
            coordinates.source_coordinate_dim,
            Some(CoordinateDimension::ThreeD)
        );
        let second = edit
            .add_default_atom_with_source_state(false, &mut bookmarks, &mut coordinates)
            .unwrap();
        assert_eq!(second, AtomId::new(2));
        assert_eq!(bookmarks.get(&-0xBADBEEF), Some(&vec![atom]));
        assert_eq!(coordinates.conformers_2d[0].coordinates()[2], [0.0; 2]);
        assert_eq!(coordinates.conformers_3d[0].coordinates()[2], [0.0; 3]);
        let (result, mapping) = edit.finish().unwrap();
        assert_eq!(result.atoms.len(), 3);
        assert_eq!(mapping.atoms().new_to_old().len(), 3);
        assert_eq!(source.atoms.len(), 1);
    }

    #[test]
    fn empty_graph_default_insertion_with_no_conformers_keeps_other_bookmarks() {
        let source = TopologyBlock::default();
        let mut edit = source.begin_batch_edit().unwrap();
        let mut bookmarks = BTreeMap::new();
        let mut coordinates = CoordinateBlock::default();
        assert_eq!(
            edit.add_default_atom_with_source_state(false, &mut bookmarks, &mut coordinates)
                .unwrap(),
            AtomId::new(0)
        );
        assert!(bookmarks.is_empty());
        assert!(coordinates.conformers_2d.is_empty());
        assert!(coordinates.conformers_3d.is_empty());
        assert_eq!(edit.working.atoms[0].element(), Element::DUMMY);
    }
}

#[cfg(test)]
mod source_replace_bond_complete_tests {
    use super::*;
    use crate::{
        BondOrder, Element, PropertyValue, SGroupCState, SubstanceGroupId, SubstanceGroupKind,
    };

    fn edit(order: BondOrder) -> TopologyBatchEdit {
        let atoms = (0..4)
            .map(|i| {
                Atom::from_spec(
                    AtomId::new(i),
                    AtomSpec::new(Element::C).with_explicit_hydrogens(if i == 0 { 3 } else { 1 }),
                )
            })
            .collect();
        let bonds = (0..3)
            .map(|i| {
                Bond::from_spec(
                    BondId::new(i),
                    BondSpec::new(
                        AtomId::new(i),
                        AtomId::new(i + 1),
                        if i == 0 { order } else { BondOrder::Single },
                    ),
                )
            })
            .collect();
        TopologyBlock::try_from_parts(atoms, bonds, vec![], vec![])
            .unwrap()
            .into_batch_edit()
            .unwrap()
    }
    fn replacement(order: BondOrder) -> Bond {
        Bond::from_spec(
            BondId::new(99),
            BondSpec::new(AtomId::new(3), AtomId::new(2), order),
        )
    }
    fn number(order: BondOrder) -> Result<f64, TopologyEditError> {
        // Inject the relevant known values at this MODEL boundary; production
        // consumers must supply the existing CORE getBondTypeAsDouble reader.
        Ok(match order {
            BondOrder::Single => 1.0,
            BondOrder::Double => 2.0,
            BondOrder::Triple => 3.0,
            BondOrder::OneAndHalf => 1.5,
            _ => panic!("unexpected fixture order"),
        })
    }
    fn uint(value: &PropertyValue) -> Result<u32, TopologyEditError> {
        value
            .as_uint()
            .map_err(|_| TopologyEditError::SourceGroupParentCycle {
                group: SubstanceGroupId::new(999),
            })
    }
    fn group(id: usize) -> SubstanceGroup {
        SubstanceGroup::new(SubstanceGroupId::new(id), SubstanceGroupKind::Data)
    }
    #[test]
    fn source_replace_bond_fractional_hydrogens_endpoints_and_dictionary_replacement() {
        let mut edit = edit(BondOrder::Single);
        edit.bond_mut(BondId::new(0))
            .unwrap()
            .set_prop("old", PropertyValue::UInt(7))
            .unwrap();
        edit.bond_mut(BondId::new(0))
            .unwrap()
            .set_computed_prop("computed", PropertyValue::Bool(true))
            .unwrap();
        let old = edit.working.bonds[0].clone();
        let mut incoming = replacement(BondOrder::OneAndHalf);
        incoming
            .set_prop("incoming", PropertyValue::UInt(9))
            .unwrap();
        incoming.set_temporary_flags(123);
        let before = incoming.clone();
        let mut reads = vec![];
        edit.replace_bond_source(
            BondId::new(0),
            &incoming,
            true,
            true,
            |order| {
                reads.push(order);
                number(order)
            },
            uint,
        )
        .unwrap();
        assert_eq!(reads, [BondOrder::OneAndHalf, BondOrder::Single]);
        assert_eq!(edit.working.atoms[0].explicit_hydrogens(), 2);
        assert_eq!(edit.working.atoms[1].explicit_hydrogens(), 0);
        let result = &edit.working.bonds[0];
        assert_eq!(
            (result.id(), result.begin(), result.end()),
            (BondId::new(0), AtomId::new(0), AtomId::new(1))
        );
        assert_eq!(result.props(), old.props());
        assert_eq!(
            result.computed_prop_names().unwrap(),
            old.computed_prop_names().unwrap()
        );
        assert_eq!(result.temporary_flags(), 123);
        assert_eq!(incoming, before);
        assert_eq!(
            edit.finish().unwrap().1.bonds.old_to_new[0],
            Some(BondId::new(0))
        );
    }
    #[test]
    fn source_replace_bond_nonincrease_and_no_preserve_keep_source_input_properties() {
        for new_order in [BondOrder::Single, BondOrder::Double] {
            let mut edit = edit(BondOrder::Double);
            let mut incoming = replacement(new_order);
            incoming
                .set_prop("incoming", PropertyValue::UInt(9))
                .unwrap();
            edit.replace_bond_source(BondId::new(0), &incoming, false, true, number, uint)
                .unwrap();
            assert_eq!(edit.working.atoms[0].explicit_hydrogens(), 3);
            assert_eq!(edit.working.atoms[1].explicit_hydrogens(), 1);
            assert_eq!(edit.working.bonds[0].props(), incoming.props());
        }
        let mut edit = edit(BondOrder::Single);
        edit.working.atoms[0].set_explicit_hydrogens(0);
        edit.replace_bond_source(
            BondId::new(0),
            &replacement(BondOrder::Triple),
            false,
            true,
            number,
            uint,
        )
        .unwrap();
        assert_eq!(edit.working.atoms[0].explicit_hydrogens(), 0);
        assert_eq!(edit.working.atoms[1].explicit_hydrogens(), 0);
    }
    #[test]
    fn source_replace_bond_group_descendants_removed_and_retained_bond_indices_decrement() {
        let mut edit = edit(BondOrder::Single);
        let mut parent = group(0).with_bonds(vec![BondId::new(0)]);
        parent.set_prop("index", PropertyValue::UInt(10)).unwrap();
        let mut child = group(1);
        child.set_prop("index", PropertyValue::UInt(20)).unwrap();
        child.set_prop("PARENT", PropertyValue::UInt(10)).unwrap();
        let mut grandchild = group(2);
        grandchild
            .set_prop("PARENT", PropertyValue::UInt(20))
            .unwrap();
        let mut orphan = group(3)
            .with_bonds(vec![BondId::new(2), BondId::new(1), BondId::new(2)])
            .with_cstates(vec![SGroupCState::new(BondId::new(2), [-0.0, 2.0, 3.0])]);
        orphan.set_prop("PARENT", PropertyValue::UInt(999)).unwrap();
        orphan
            .set_prop("raw", PropertyValue::String("x".into()))
            .unwrap();
        edit.working.substance_groups = vec![parent, child, grandchild, orphan];
        edit.replace_bond_source(
            BondId::new(0),
            &replacement(BondOrder::Double),
            false,
            false,
            number,
            uint,
        )
        .unwrap();
        assert_eq!(edit.working.substance_groups.len(), 1);
        let surviving = &edit.working.substance_groups[0];
        assert_eq!(surviving.id(), SubstanceGroupId::new(0));
        assert_eq!(
            surviving.bonds(),
            [BondId::new(1), BondId::new(0), BondId::new(1)]
        );
        assert_eq!(surviving.cstates()[0].bond(), BondId::new(1));
        assert_eq!(
            surviving.props().get(b"PARENT".as_slice()),
            Some(&PropertyValue::UInt(999))
        );
    }
    #[test]
    fn source_replace_bond_no_removal_skips_index_conversion_keep_groups_skips_all_group_reads() {
        let mut edit = edit(BondOrder::Single);
        let mut retained = group(0).with_bonds(vec![BondId::new(2)]);
        retained
            .set_prop("index", PropertyValue::Bool(false))
            .unwrap();
        retained.set_prop("PARENT", PropertyValue::UInt(2)).unwrap();
        edit.working.substance_groups.push(retained);
        edit.replace_bond_source(
            BondId::new(0),
            &replacement(BondOrder::Double),
            false,
            false,
            number,
            uint,
        )
        .unwrap();
        assert_eq!(edit.working.substance_groups[0].bonds(), [BondId::new(1)]);
        let before = edit.working.substance_groups.clone();
        edit.replace_bond_source(
            BondId::new(0),
            &replacement(BondOrder::Double),
            false,
            true,
            number,
            |_| panic!("keepSGroups must not read group properties"),
        )
        .unwrap();
        assert_eq!(edit.working.substance_groups, before);
    }
    #[test]
    fn source_replace_bond_range_order_reader_and_group_conversion_failures_retain_source_effects()
    {
        let mut edit = edit(BondOrder::Single);
        let before = edit.working.clone();
        let error = edit
            .replace_bond_source::<TopologyEditError>(
                BondId::new(9),
                &replacement(BondOrder::Double),
                false,
                false,
                |_| panic!("range before type read"),
                |_| panic!("range before group read"),
            )
            .unwrap_err();
        assert!(matches!(error, TopologyEditError::BondOutOfRange { .. }));
        assert_eq!(edit.working, before);
        let mut read_count = 0;
        assert!(
            edit.replace_bond_source(
                BondId::new(0),
                &replacement(BondOrder::Double),
                false,
                false,
                |_| {
                    read_count += 1;
                    Err(TopologyEditError::SourceGroupParentCycle {
                        group: SubstanceGroupId::new(0),
                    })
                },
                uint
            )
            .is_err()
        );
        assert_eq!(read_count, 1);
        assert_eq!(edit.working, before);
        let mut bad = group(1);
        bad.set_prop("PARENT", PropertyValue::Bool(false)).unwrap();
        edit.working.substance_groups = vec![group(0), bad];
        assert!(
            edit.replace_bond_source(
                BondId::new(0),
                &replacement(BondOrder::Double),
                false,
                false,
                number,
                uint
            )
            .is_err()
        );
        assert_eq!(edit.working.bonds[0].order(), BondOrder::Double);
        assert_eq!(edit.working.atoms[0].explicit_hydrogens(), 2);
        assert_eq!(edit.working.atoms[1].explicit_hydrogens(), 0);
        assert_eq!(edit.working.substance_groups.len(), 2);
    }
}

/// Explicit transient detached source state; neither branch owns live state.
#[allow(dead_code)]
enum SourceBondRemovalState<'a, E> {
    Batch {
        pending: &'a mut Vec<bool>,
    },
    Immediate {
        bookmarks: &'a mut std::collections::BTreeMap<i32, Vec<BondId>>,
        reset_ring: &'a mut dyn FnMut(),
        uint_reader: &'a mut dyn FnMut(&crate::PropertyValue) -> Result<u32, E>,
    },
}

#[allow(dead_code)]
fn remove_bond_source<E: From<TopologyEditError>>(
    topology: &mut TopologyBlock,
    begin: AtomId,
    end: AtomId,
    state: SourceBondRemovalState<'_, E>,
) -> Result<(), E> {
    // RDKit❗❌: void RWMol::removeBond(unsigned int aid1, unsigned int aid2) {
    // RDKit❗❌:   URANGE_CHECK(aid1, getNumAtoms());
    // RDKit❗❌:   URANGE_CHECK(aid2, getNumAtoms());
    // RDKit❗❌:   auto *bnd = getBondBetweenAtoms(aid1, aid2);
    // RDKit❗❌:   if (!bnd) {
    // RDKit❗❌:     return;
    // RDKit❗❌:   }
    // RDKit❗❌:   auto idx = bnd->getIdx();
    // RDKit❗❌:   if (dp_delBonds) {
    // RDKit❗❌:     // we're in a batch edit
    // RDKit❗❌:     // if bonds have been added since we started, resize dp_delBonds
    // RDKit❗❌:     if (dp_delBonds->size() < getNumBonds()) {
    // RDKit❗❌:       dp_delBonds->resize(getNumBonds());
    // RDKit❗❌:     }
    // RDKit❗❌:     dp_delBonds->set(idx);
    // RDKit❗❌:     return;
    // RDKit❗❌:   }
    // RDKit❗❌:
    // RDKit❗❌:   // remove any bookmarks which point to this bond:
    // RDKit❗❌:   BOND_BOOKMARK_MAP *marks = getBondBookmarks();
    // RDKit❗❌:   auto markI = marks->begin();
    // RDKit❗❌:   while (markI != marks->end()) {
    // RDKit❗❌:     BOND_PTR_LIST &bonds = markI->second;
    // RDKit❗❌:     // we need to copy the iterator then increment it, because the
    // RDKit❗❌:     // deletion we're going to do in clearBondBookmark will invalidate
    // RDKit❗❌:     // it.
    // RDKit❗❌:     auto tmpI = markI;
    // RDKit❗❌:     ++markI;
    // RDKit❗❌:     if (std::find(bonds.begin(), bonds.end(), bnd) != bonds.end()) {
    // RDKit❗❌:       clearBondBookmark(tmpI->first, bnd);
    // RDKit❗❌:     }
    // RDKit❗❌:   }
    // RDKit❗❌:
    // RDKit❗❌:   // loop over neighboring double bonds and remove their stereo atom
    // RDKit❗❌:   //  information. This is definitely now invalid (was github issue 8)
    // RDKit❗❌:   auto beginAtm = bnd->getBeginAtom();
    // RDKit❗❌:   auto endAtm = bnd->getEndAtom();
    // RDKit❗❌:   std::vector<std::vector<Atom *>> bond_atoms = {{beginAtm, endAtm},
    // RDKit❗❌:                                                  {endAtm, beginAtm}};
    // RDKit❗❌:   for (const auto &atoms : bond_atoms) {
    // RDKit❗❌:     for (auto obnd : this->atomBonds(atoms[0])) {
    // RDKit❗❌:       if (obnd == bnd) {
    // RDKit❗❌:         continue;
    // RDKit❗❌:       }
    // RDKit❗❌:       if (std::find(obnd->getStereoAtoms().begin(),
    // RDKit❗❌:                     obnd->getStereoAtoms().end(),
    // RDKit❗❌:                     atoms[1]->getIdx()) != obnd->getStereoAtoms().end()) {
    // RDKit❗❌:         // github #6900 if we remove stereo atoms we need to remove
    // RDKit❗❌:         //  the CIS and or TRANS since this requires stereo atoms
    // RDKit❗❌:         if (obnd->getStereo() == Bond::BondStereo::STEREOCIS ||
    // RDKit❗❌:             obnd->getStereo() == Bond::BondStereo::STEREOTRANS) {
    // RDKit❗❌:           obnd->setStereo(Bond::BondStereo::STEREONONE);
    // RDKit❗❌:         }
    // RDKit❗❌:
    // RDKit❗❌:         obnd->getStereoAtoms().clear();
    // RDKit❗❌:       }
    // RDKit❗❌:     }
    // RDKit❗❌:   }
    // RDKit❗❌:   // reset our ring info structure, because it is pretty likely
    // RDKit❗❌:   // to be wrong now:
    // RDKit❗❌:   dp_ringInfo->reset();
    // RDKit❗❌:
    // RDKit❗❌:   removeSubstanceGroupsReferencingBond(*this, idx);
    // RDKit❗❌:   removeBondFromGroups(bnd, d_stereo_groups);
    // RDKit❗❌:
    // RDKit❗❌:   // loop over all bonds with higher indices and update their indices
    // RDKit❗❌:   for (auto bond : bonds()) {
    // RDKit❗❌:     if (bond->getIdx() > idx) {
    // RDKit❗❌:       bond->setIdx(bond->getIdx() - 1);
    // RDKit❗❌:     }
    // RDKit❗❌:   }
    // RDKit❗❌:
    // RDKit❗❌:   auto vd1 = boost::vertex(bnd->getBeginAtomIdx(), d_graph);
    // RDKit❗❌:   auto vd2 = boost::vertex(bnd->getEndAtomIdx(), d_graph);
    // RDKit❗❌:   boost::remove_edge(vd1, vd2, d_graph);
    // RDKit❗❌:   delete bnd;
    // RDKit❗❌:   --numBonds;
    // RDKit❗❌: }
    // Detached IDs cannot model dangling source bookmark/stereo pointers if
    // duplicate references survive the native first-only erasure. Such native
    // invalid pointer states remain a behavior gap, not silently cleared.
    // Cost: canonical row lookup/incident scans and final CSR rebuild add O(E)
    // work versus source adjacency-local lookup. No whole block/property clone.
    for atom in [begin, end] {
        if atom.index() >= topology.atoms.len() {
            return Err(TopologyEditError::AtomOutOfRange {
                atom,
                atom_count: topology.atoms.len(),
            }
            .into());
        }
    }
    let Some(index) = topology.bonds.iter().position(|b| {
        (b.begin() == begin && b.end() == end) || (b.begin() == end && b.end() == begin)
    }) else {
        return Ok(());
    };
    let bond_id = topology.bonds[index].id();
    let (bookmarks, reset_ring, uint_reader) = match state {
        SourceBondRemovalState::Batch { pending } => {
            if pending.len() < topology.bonds.len() {
                pending.resize(topology.bonds.len(), false);
            }
            pending[index] = true;
            return Ok(());
        }
        SourceBondRemovalState::Immediate {
            bookmarks,
            reset_ring,
            uint_reader,
        } => (bookmarks, reset_ring, uint_reader),
    };
    // Native clearBondBookmark removes the first matching entry, then removes
    // an empty map key. It deliberately does not erase all repeated pointers.
    clear_source_bond_bookmarks(bookmarks, bond_id);
    clear_source_incident_bond_stereo(topology, index).map_err(E::from)?;
    reset_ring();
    crate::sgroup::remove_source_groups_referencing_bond(
        &mut topology.substance_groups,
        bond_id,
        uint_reader,
    )?;
    for group in &mut topology.stereo_groups {
        group.remove_bond(bond_id);
    }
    topology.stereo_groups.retain(|group| !group.is_empty());
    for bond in &mut topology.bonds {
        if bond.id().index() > index {
            bond.set_id_for_construction(BondId::new(bond.id().index() - 1));
        }
    }
    // Pointer references track their native object across index changes.
    for values in bookmarks.values_mut() {
        for id in values {
            if id.index() > index {
                *id = BondId::new(id.index() - 1);
            }
        }
    }
    for group in &mut topology.stereo_groups {
        group.shift_source_bond_ids_after_removed(bond_id);
    }
    topology.bonds.remove(index);
    topology.adjacency = AdjacencyList::from_topology(topology.atoms.len(), &topology.bonds);
    Ok(())
}

#[cfg(test)]
mod source_remove_bond_complete_tests {
    use super::*;
    use crate::{
        BondOrder, Element, PropertyValue, StereoGroupKind, SubstanceGroupId, SubstanceGroupKind,
    };
    use std::{cell::Cell, collections::BTreeMap};
    fn topology(stereo: BondStereo) -> TopologyBlock {
        let atoms = (0..4)
            .map(|i| Atom::from_spec(AtomId::new(i), AtomSpec::new(Element::C)))
            .collect();
        let mut middle = Bond::from_spec(
            BondId::new(1),
            BondSpec::new(AtomId::new(1), AtomId::new(2), BondOrder::Double),
        );
        middle.set_stereo_atoms(Some([AtomId::new(0), AtomId::new(3)]));
        middle.set_stereo(stereo).unwrap();
        TopologyBlock::try_from_parts(
            atoms,
            vec![
                Bond::from_spec(
                    BondId::new(0),
                    BondSpec::new(AtomId::new(0), AtomId::new(1), BondOrder::Single),
                ),
                middle,
                Bond::from_spec(
                    BondId::new(2),
                    BondSpec::new(AtomId::new(2), AtomId::new(3), BondOrder::Single),
                ),
            ],
            vec![],
            vec![],
        )
        .unwrap()
    }
    fn uint(value: &PropertyValue) -> Result<u32, TopologyEditError> {
        value
            .as_uint()
            .map_err(|_| TopologyEditError::SourceGroupParentCycle {
                group: SubstanceGroupId::new(999),
            })
    }
    #[test]
    fn source_remove_bond_batch_mask_resizes_and_retains_all_other_state() {
        let mut topology = topology(BondStereo::Cis);
        let before = topology.clone();
        let mut pending = vec![];
        remove_bond_source::<TopologyEditError>(
            &mut topology,
            AtomId::new(1),
            AtomId::new(0),
            SourceBondRemovalState::Batch {
                pending: &mut pending,
            },
        )
        .unwrap();
        assert_eq!(pending, [true, false, false]);
        assert_eq!(topology, before);
        remove_bond_source::<TopologyEditError>(
            &mut topology,
            AtomId::new(0),
            AtomId::new(1),
            SourceBondRemovalState::Batch {
                pending: &mut pending,
            },
        )
        .unwrap();
        assert_eq!(pending, [true, false, false]);
        let mut edit = topology.into_batch_edit().unwrap();
        let new_atom = edit.add_atom(AtomSpec::new(Element::N));
        edit.add_bond(BondSpec::new(new_atom, AtomId::new(3), BondOrder::Single))
            .unwrap();
        edit.remove_bond_between_atoms(new_atom, AtomId::new(3))
            .unwrap();
        assert!(edit.remove_bonds[3]);
        assert_eq!(edit.working.bonds.len(), 4);
    }
    #[test]
    fn source_remove_bond_checks_endpoint_order_and_absent_edge_is_noop() {
        let mut topology = topology(BondStereo::Cis);
        let before = topology.clone();
        let mut pending = vec![];
        remove_bond_source::<TopologyEditError>(
            &mut topology,
            AtomId::new(0),
            AtomId::new(3),
            SourceBondRemovalState::Batch {
                pending: &mut pending,
            },
        )
        .unwrap();
        assert!(pending.is_empty());
        let error = remove_bond_source::<TopologyEditError>(
            &mut topology,
            AtomId::new(9),
            AtomId::new(8),
            SourceBondRemovalState::Batch {
                pending: &mut pending,
            },
        )
        .unwrap_err();
        assert!(
            matches!(error, TopologyEditError::AtomOutOfRange { atom, .. } if atom == AtomId::new(9))
        );
        assert_eq!(topology, before);
        assert!(pending.is_empty());
    }
    #[test]
    fn source_remove_bond_immediate_stereo_bookmark_groups_and_dense_indices() {
        for stereo in [
            BondStereo::Cis,
            BondStereo::Trans,
            BondStereo::E,
            BondStereo::Z,
            BondStereo::Any,
        ] {
            let mut topology = topology(stereo);
            topology.stereo_groups = vec![
                StereoGroup::new(StereoGroupKind::Absolute, vec![], vec![BondId::new(0)])
                    .expect("valid distinct stereo members"),
                StereoGroup::new(
                    StereoGroupKind::And,
                    vec![AtomId::new(1)],
                    vec![BondId::new(2)],
                )
                .expect("valid distinct stereo members"),
            ];
            topology.substance_groups = vec![
                SubstanceGroup::new(SubstanceGroupId::new(0), SubstanceGroupKind::Data)
                    .with_bonds(vec![BondId::new(0)]),
                SubstanceGroup::new(SubstanceGroupId::new(1), SubstanceGroupKind::Data)
                    .with_bonds(vec![BondId::new(2)]),
            ];
            let mut bookmarks = BTreeMap::from([
                (1, vec![BondId::new(0)]),
                (2, vec![BondId::new(2), BondId::new(0)]),
            ]);
            let resets = Cell::new(0);
            remove_bond_source(
                &mut topology,
                AtomId::new(1),
                AtomId::new(0),
                SourceBondRemovalState::Immediate {
                    bookmarks: &mut bookmarks,
                    reset_ring: &mut || resets.set(resets.get() + 1),
                    uint_reader: &mut uint,
                },
            )
            .unwrap();
            assert_eq!(resets.get(), 1);
            assert_eq!(bookmarks, BTreeMap::from([(2, vec![BondId::new(1)])]));
            assert_eq!(topology.bonds.len(), 2);
            assert_eq!(topology.bonds[0].id(), BondId::new(0));
            assert_eq!(topology.bonds[0].stereo_atoms(), None);
            assert_eq!(
                topology.bonds[0].stereo(),
                if matches!(stereo, BondStereo::Cis | BondStereo::Trans) {
                    BondStereo::None
                } else {
                    stereo
                }
            );
            assert_eq!(topology.substance_groups[0].bonds(), [BondId::new(1)]);
            assert_eq!(topology.stereo_groups.len(), 1);
            assert_eq!(topology.stereo_groups[0].bonds(), [BondId::new(1)]);
            topology.validate().unwrap();
        }
    }
    #[test]
    fn source_remove_bond_group_error_occurs_after_bookmark_stereo_and_ring_effects() {
        let mut topology = topology(BondStereo::Cis);
        let mut bad = SubstanceGroup::new(SubstanceGroupId::new(0), SubstanceGroupKind::Data);
        bad.set_prop("PARENT", PropertyValue::Bool(false)).unwrap();
        topology.substance_groups.push(bad);
        let mut bookmarks = BTreeMap::from([(1, vec![BondId::new(0), BondId::new(0)])]);
        let resets = Cell::new(0);
        let error = remove_bond_source(
            &mut topology,
            AtomId::new(0),
            AtomId::new(1),
            SourceBondRemovalState::Immediate {
                bookmarks: &mut bookmarks,
                reset_ring: &mut || resets.set(1),
                uint_reader: &mut |value| {
                    assert_eq!(resets.get(), 1);
                    uint(value)
                },
            },
        )
        .unwrap_err();
        assert!(matches!(
            error,
            TopologyEditError::SourceGroupParentCycle { .. }
        ));
        assert_eq!(bookmarks[&1], [BondId::new(0)]); // native first-only erasure
        assert_eq!(topology.bonds[1].stereo_atoms(), None);
        assert_eq!(topology.bonds[1].stereo(), BondStereo::None);
        assert_eq!(topology.bonds.len(), 3); // graph deletion has not been reached
        assert_eq!(topology.bonds[2].id(), BondId::new(2));
    }
}

fn update_source_bond_endpoints_property(
    bond: &mut Bond,
    removed_atom_index: usize,
    source_property: &crate::PropertyText,
) -> Result<(), TopologyEditError> {
    let bond_id = bond.id();
    // RDKit❗❌:       if ('(' == sprop.front() && ')' == sprop.back()) {
    // RDKit❗❌:         sprop = sprop.substr(1, sprop.length() - 2);
    // RDKit❗❌:
    // RDKit❗❌:         // This is doing what ParseV3000Array would do.
    // RDKit❗❌:         boost::char_separator<char> sep(" ");
    // RDKit❗❌:         boost::tokenizer<boost::char_separator<char>> tokens(sprop, sep);
    // RDKit❗❌:         unsigned int num_ats = std::stod(*tokens.begin());
    // RDKit❗❌:         std::vector<unsigned int> oats;
    // RDKit❗❌:         auto beg = tokens.begin();
    // RDKit❗❌:         ++beg;
    // RDKit❗❌:         std::transform(beg, tokens.end(), std::back_inserter(oats),
    // RDKit❗❌:                        [](const std::string &a) { return std::stod(a); });
    // RDKit❗❌:
    // RDKit❗❌:         auto idx_pos = std::find(oats.begin(), oats.end(), idx + 1);
    // RDKit❗❌:         if (idx_pos != oats.end()) {
    // RDKit❗❌:           oats.erase(idx_pos);
    // RDKit❗❌:           --num_ats;
    // RDKit❗❌:         }
    // RDKit❗❌:         if (!num_ats) {
    // RDKit❗❌:           bond->clearProp(RDKit::common_properties::_MolFileBondEndPts);
    // RDKit❗❌:           bond->clearProp(common_properties::_MolFileBondAttach);
    // RDKit❗❌:         } else {
    // RDKit❗❌:           sprop = "(" + std::to_string(num_ats) + " ";
    // RDKit❗❌:           for (auto &i : oats) {
    // RDKit❗❌:             if (i > idx + 1) {
    // RDKit❗❌:               --i;
    // RDKit❗❌:             }
    // RDKit❗❌:             sprop += std::to_string(i) + " ";
    // RDKit❗❌:           }
    // RDKit❗❌:           sprop[sprop.length() - 1] = ')';
    // RDKit❗❌:           bond->setProp(RDKit::common_properties::_MolFileBondEndPts, sprop);
    // RDKit❗❌:         }
    // RDKit❗❌:       }
    // RDKit❗❌:     }
    let Some(contents) = source_property
        .as_bytes()
        .strip_prefix(b"(")
        .and_then(|value| value.strip_suffix(b")"))
    else {
        return Ok(());
    };

    let mut tokens = contents
        .split(|&byte| byte == b' ')
        .filter(|token| !token.is_empty());
    let Some(count_token) = tokens.next() else {
        return Err(TopologyEditError::BondEndPointsParse {
            bond: bond_id,
            token_index: 0,
            token: crate::PropertyText::new(),
            kind: BondEndPointsParseErrorKind::InvalidArgument,
        });
    };
    let parse_token = |token: &[u8], token_index: usize| -> Result<u32, TopologyEditError> {
        let (value, _) = parse_endpts_decimal_prefix(token).map_err(|kind| {
            TopologyEditError::BondEndPointsParse {
                bond: bond_id,
                token_index,
                token: token.into(),
                kind,
            }
        })?;
        endpts_value_to_u32(value).map_err(|kind| TopologyEditError::BondEndPointsParse {
            bond: bond_id,
            token_index,
            token: token.into(),
            kind,
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
        bond.clear_prop("_MolFileBondEndPts")
            .map_err(TopologyEditError::InvalidBond)?;
        bond.clear_prop("_MolFileBondAttach")
            .map_err(TopologyEditError::InvalidBond)?;
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
    Ok(())
}

#[allow(dead_code)]
enum SourceAtomRemovalState<'a, E> {
    Batch {
        atoms: &'a mut Vec<bool>,
        bonds: &'a mut Vec<bool>,
    },
    Immediate {
        atom_bookmarks: &'a mut std::collections::BTreeMap<i32, Vec<AtomId>>,
        bond_bookmarks: &'a mut std::collections::BTreeMap<i32, Vec<BondId>>,
        coordinates: &'a mut crate::CoordinateBlock,
        properties: &'a mut crate::MoleculeProperties,
        reset_ring: &'a mut dyn FnMut(),
        uint_reader: &'a mut dyn FnMut(&crate::PropertyValue) -> Result<u32, E>,
        text_reader: &'a mut dyn FnMut(&crate::PropertyValue) -> Result<crate::PropertyText, E>,
    },
}

#[allow(dead_code)]
fn remove_atom_source_index<E: From<TopologyEditError>>(
    topology: &mut TopologyBlock,
    index: AtomId,
    state: SourceAtomRemovalState<'_, E>,
) -> Result<(), E> {
    // RDKit❗✔️: void RWMol::removeAtom(unsigned int idx) { removeAtom(getAtomWithIdx(idx)); }
    // RDKit❗✔️: Atom *ROMol::getAtomWithIdx(unsigned int idx) {
    // RDKit❗✔️:   URANGE_CHECK(idx, getNumAtoms());
    // RDKit❗✔️:   auto vd = boost::vertex(idx, d_graph);
    // RDKit❗✔️:   auto res = d_graph[vd];
    // RDKit❗✔️:   POSTCONDITION(res, "");
    // RDKit❗✔️:   return res;
    // RDKit❗✔️: }
    // The indexed row supplies the source atom identity. No copy or traversal;
    // source uint32/null/ownership and dependency gaps remain explicit.
    let atom = topology
        .atoms
        .get(index.index())
        .ok_or_else(|| {
            E::from(TopologyEditError::AtomOutOfRange {
                atom: index,
                atom_count: topology.atoms.len(),
            })
        })?
        .id();
    remove_atom_source_default(topology, atom, state)
}

#[allow(dead_code)]
fn remove_atom_source_default<E: From<TopologyEditError>>(
    topology: &mut TopologyBlock,
    atom: AtomId,
    state: SourceAtomRemovalState<'_, E>,
) -> Result<(), E> {
    // RDKit❗✔️: void RWMol::removeAtom(Atom *atom) { removeAtom(atom, true); }
    // The source literal true reaches the sole complete bool overload. Native
    // pointer/ownership and source-state gaps are retained at that dependency.
    // Cost: constant forwarding, with no allocation, copy or extra traversal.
    remove_atom_source(topology, atom, true, state)
}

#[allow(dead_code)]
fn remove_atom_source<E: From<TopologyEditError>>(
    topology: &mut TopologyBlock,
    atom: AtomId,
    clear_props: bool,
    state: SourceAtomRemovalState<'_, E>,
) -> Result<(), E> {
    // RDKit❗❌: void RWMol::removeAtom(Atom *atom, bool clearProps) {
    // RDKit❗❌:   PRECONDITION(atom, "NULL atom provided");
    // RDKit❗❌:   PRECONDITION(static_cast<RWMol *>(&atom->getOwningMol()) == this,
    // RDKit❗❌:                "atom not owned by this molecule");
    // RDKit❗❌:   unsigned int idx = atom->getIdx();
    // RDKit❗❌:
    // RDKit❗❌:   // remove bonds attached to the atom
    // RDKit❗❌:   //  In batch mode this will schedule bond removal
    // RDKit❗❌:   std::vector<std::pair<unsigned int, unsigned int>> nbrs;
    // RDKit❗❌:   ADJ_ITER b1, b2;
    // RDKit❗❌:   boost::tie(b1, b2) = getAtomNeighbors(atom);
    // RDKit❗❌:   while (b1 != b2) {
    // RDKit❗❌:     nbrs.emplace_back(atom->getIdx(), rdcast<unsigned int>(*b1));
    // RDKit❗❌:     ++b1;
    // RDKit❗❌:   }
    // RDKit❗❌:   for (auto &nbr : nbrs) {
    // RDKit❗❌:     removeBond(nbr.first, nbr.second);
    // RDKit❗❌:   }
    // RDKit❗❌:
    // RDKit❗❌:   if (dp_delAtoms) {
    // RDKit❗❌:     // we're in a batch edit
    // RDKit❗❌:     // if atoms have been added since we started, resize dp_delAtoms
    // RDKit❗❌:     if (dp_delAtoms->size() < getNumAtoms()) {
    // RDKit❗❌:       dp_delAtoms->resize(getNumAtoms());
    // RDKit❗❌:     }
    // RDKit❗❌:     dp_delAtoms->set(idx);
    // RDKit❗❌:     return;
    // RDKit❗❌:   }
    // RDKit❗❌:
    // RDKit❗❌:   // remove any bookmarks which point to this atom:
    // RDKit❗❌:   ATOM_BOOKMARK_MAP *marks = getAtomBookmarks();
    // RDKit❗❌:   auto markI = marks->begin();
    // RDKit❗❌:   while (markI != marks->end()) {
    // RDKit❗❌:     const ATOM_PTR_LIST &atoms = markI->second;
    // RDKit❗❌:     // we need to copy the iterator then increment it, because the
    // RDKit❗❌:     // deletion we're going to do in clearAtomBookmark will invalidate
    // RDKit❗❌:     // it.
    // RDKit❗❌:     auto tmpI = markI;
    // RDKit❗❌:     ++markI;
    // RDKit❗❌:     if (std::find(atoms.begin(), atoms.end(), atom) != atoms.end()) {
    // RDKit❗❌:       clearAtomBookmark(tmpI->first, atom);
    // RDKit❗❌:     }
    // RDKit❗❌:   }
    // RDKit❗❌:
    // RDKit❗❌:   // loop over all atoms with higher indices and update their indices
    // RDKit❗❌:   for (unsigned int i = idx + 1; i < getNumAtoms(); i++) {
    // RDKit❗❌:     Atom *higher_index_atom = getAtomWithIdx(i);
    // RDKit❗❌:     higher_index_atom->setIdx(i - 1);
    // RDKit❗❌:   }
    // RDKit❗❌:
    // RDKit❗❌:   // do the same with the coordinates in the conformations
    // RDKit❗❌:   for (auto conf : d_confs) {
    // RDKit❗❌:     RDGeom::POINT3D_VECT &positions = conf->getPositions();
    // RDKit❗❌:     auto pi = positions.begin();
    // RDKit❗❌:     for (unsigned int i = 0; i < getNumAtoms() - 1; i++) {
    // RDKit❗❌:       ++pi;
    // RDKit❗❌:       if (i >= idx) {
    // RDKit❗❌:         positions[i] = positions[i + 1];
    // RDKit❗❌:       }
    // RDKit❗❌:     }
    // RDKit❗❌:     positions.erase(pi);
    // RDKit❗❌:   }
    // RDKit❗❌:   // now deal with bonds:
    // RDKit❗❌:   //   their end indices may need to be decremented and their
    // RDKit❗❌:   //   indices will need to be handled and if they have an
    // RDKit❗❌:   //   ENDPTS prop that includes idx, it will need updating.
    // RDKit❗❌:   unsigned int nBonds = 0;
    // RDKit❗❌:   EDGE_ITER beg, end;
    // RDKit❗❌:   boost::tie(beg, end) = getEdges();
    // RDKit❗❌:   std::string sprop;
    // RDKit❗❌:   while (beg != end) {
    // RDKit❗❌:     Bond *bond = d_graph[*beg++];
    // RDKit❗❌:     if (bond->getPropIfPresent(RDKit::common_properties::_MolFileBondEndPts,
    // RDKit❗❌:                                sprop)) {
    // RDKit❗❌:       // This would ideally use ParseV3000Array but I'm buggered if I can get
    // RDKit❗❌:       // the linker to find it.
    // RDKit❗❌:       //      std::vector<unsigned int> oats =
    // RDKit❗❌:       //          RDKit::SGroupParsing::ParseV3000Array<unsigned int>(sprop);
    // RDKit❗❌:       if ('(' == sprop.front() && ')' == sprop.back()) {
    // RDKit❗❌:         sprop = sprop.substr(1, sprop.length() - 2);
    // RDKit❗❌:
    // RDKit❗❌:         // This is doing what ParseV3000Array would do.
    // RDKit❗❌:         boost::char_separator<char> sep(" ");
    // RDKit❗❌:         boost::tokenizer<boost::char_separator<char>> tokens(sprop, sep);
    // RDKit❗❌:         unsigned int num_ats = std::stod(*tokens.begin());
    // RDKit❗❌:         std::vector<unsigned int> oats;
    // RDKit❗❌:         auto beg = tokens.begin();
    // RDKit❗❌:         ++beg;
    // RDKit❗❌:         std::transform(beg, tokens.end(), std::back_inserter(oats),
    // RDKit❗❌:                        [](const std::string &a) { return std::stod(a); });
    // RDKit❗❌:
    // RDKit❗❌:         auto idx_pos = std::find(oats.begin(), oats.end(), idx + 1);
    // RDKit❗❌:         if (idx_pos != oats.end()) {
    // RDKit❗❌:           oats.erase(idx_pos);
    // RDKit❗❌:           --num_ats;
    // RDKit❗❌:         }
    // RDKit❗❌:         if (!num_ats) {
    // RDKit❗❌:           bond->clearProp(RDKit::common_properties::_MolFileBondEndPts);
    // RDKit❗❌:           bond->clearProp(common_properties::_MolFileBondAttach);
    // RDKit❗❌:         } else {
    // RDKit❗❌:           sprop = "(" + std::to_string(num_ats) + " ";
    // RDKit❗❌:           for (auto &i : oats) {
    // RDKit❗❌:             if (i > idx + 1) {
    // RDKit❗❌:               --i;
    // RDKit❗❌:             }
    // RDKit❗❌:             sprop += std::to_string(i) + " ";
    // RDKit❗❌:           }
    // RDKit❗❌:           sprop[sprop.length() - 1] = ')';
    // RDKit❗❌:           bond->setProp(RDKit::common_properties::_MolFileBondEndPts, sprop);
    // RDKit❗❌:         }
    // RDKit❗❌:       }
    // RDKit❗❌:     }
    // RDKit❗❌:     unsigned int tmpIdx = bond->getBeginAtomIdx();
    // RDKit❗❌:     if (tmpIdx > idx) {
    // RDKit❗❌:       bond->setBeginAtomIdx(tmpIdx - 1);
    // RDKit❗❌:     }
    // RDKit❗❌:     tmpIdx = bond->getEndAtomIdx();
    // RDKit❗❌:     if (tmpIdx > idx) {
    // RDKit❗❌:       bond->setEndAtomIdx(tmpIdx - 1);
    // RDKit❗❌:     }
    // RDKit❗❌:     bond->setIdx(nBonds++);
    // RDKit❗❌:     for (auto bsi = bond->getStereoAtoms().begin();
    // RDKit❗❌:          bsi != bond->getStereoAtoms().end(); ++bsi) {
    // RDKit❗❌:       if ((*bsi) == rdcast<int>(idx)) {
    // RDKit❗❌:         bond->getStereoAtoms().clear();
    // RDKit❗❌:         break;
    // RDKit❗❌:       } else if ((*bsi) > rdcast<int>(idx)) {
    // RDKit❗❌:         --(*bsi);
    // RDKit❗❌:       }
    // RDKit❗❌:     }
    // RDKit❗❌:   }
    // RDKit❗❌:
    // RDKit❗❌:   removeSubstanceGroupsReferencingAtom(*this, idx);
    // RDKit❗❌:
    // RDKit❗❌:   // Remove this atom from any stereo group
    // RDKit❗❌:   removeAtomFromGroups(atom, d_stereo_groups);
    // RDKit❗❌:
    // RDKit❗❌:   // clear computed properties and reset our ring info structure
    // RDKit❗❌:   // they are pretty likely to be wrong now:
    // RDKit❗❌:   if (clearProps) {
    // RDKit❗❌:     clearComputedProps(true);
    // RDKit❗❌:   }
    // RDKit❗❌:
    // RDKit❗❌:   atom->setOwningMol(nullptr);
    // RDKit❗❌:
    // RDKit❗❌:   // remove all connections to the atom:
    // RDKit❗❌:   MolGraph::vertex_descriptor vd = boost::vertex(idx, d_graph);
    // RDKit❗❌:   boost::clear_vertex(vd, d_graph);
    // RDKit❗❌:   // finally remove the vertex itself
    // RDKit❗❌:   boost::remove_vertex(vd, d_graph);
    // RDKit❗❌:   delete atom;
    // RDKit❗❌: }
    // Source owning-pointer/null predicates project to the checked canonical
    // detached row. Native invalid-pointer/width/allocator states are unmodeled.
    // Stable stereo-group IDs track native pointers on successful completion;
    // automatic pointer index changes during an earlier failure remain a gap.
    // Row scan/CSR rebuilding and copied mixed-dimension order add cost.
    let atom_count = topology.atoms.len();
    if atom.index() >= atom_count {
        return Err(TopologyEditError::AtomOutOfRange { atom, atom_count }.into());
    }
    let neighbors: Vec<_> = topology
        .bonds
        .iter()
        .filter_map(|bond| {
            if bond.begin() == atom {
                Some(bond.end())
            } else if bond.end() == atom {
                Some(bond.begin())
            } else {
                None
            }
        })
        .collect();
    let (
        atom_bookmarks,
        bond_bookmarks,
        coordinates,
        properties,
        reset_ring,
        uint_reader,
        text_reader,
    ) = match state {
        SourceAtomRemovalState::Batch { atoms, bonds } => {
            for neighbor in neighbors {
                remove_bond_source::<E>(
                    topology,
                    atom,
                    neighbor,
                    SourceBondRemovalState::Batch { pending: bonds },
                )?;
            }
            if atoms.len() < atom_count {
                atoms.resize(atom_count, false);
            }
            atoms[atom.index()] = true;
            return Ok(());
        }
        SourceAtomRemovalState::Immediate {
            atom_bookmarks,
            bond_bookmarks,
            coordinates,
            properties,
            reset_ring,
            uint_reader,
            text_reader,
        } => (
            atom_bookmarks,
            bond_bookmarks,
            coordinates,
            properties,
            reset_ring,
            uint_reader,
            text_reader,
        ),
    };
    for neighbor in neighbors {
        remove_bond_source(
            topology,
            atom,
            neighbor,
            SourceBondRemovalState::Immediate {
                bookmarks: &mut *bond_bookmarks,
                reset_ring: &mut *reset_ring,
                uint_reader: &mut *uint_reader,
            },
        )?;
    }
    clear_source_atom_bookmarks(atom_bookmarks, atom);
    for higher in topology.atoms.iter_mut().skip(atom.index() + 1) {
        higher.set_source_index(AtomId::new(higher.id().index() - 1));
    }
    for refs in atom_bookmarks.values_mut() {
        for id in refs {
            if id.index() > atom.index() {
                *id = AtomId::new(id.index() - 1);
            }
        }
    }
    // Actual source mixed-conformer append order controls failure order.
    let order = match &coordinates.source_conformer_order {
        Some(order) => order.clone(),
        None if coordinates.conformers_2d.is_empty() => {
            vec![crate::CoordinateDimension::ThreeD; coordinates.conformers_3d.len()]
        }
        None if coordinates.conformers_3d.is_empty() => {
            vec![crate::CoordinateDimension::TwoD; coordinates.conformers_2d.len()]
        }
        None => {
            return Err(TopologyEditError::Coordinate(
                crate::CoordinateValidationError::MissingSourceConformerOrder,
            )
            .into());
        }
    };
    let mut two_d = 0;
    let mut three_d = 0;
    for dimension in order {
        match dimension {
            crate::CoordinateDimension::TwoD => {
                let conformer = coordinates.conformers_2d.get_mut(two_d).ok_or_else(|| {
                    E::from(TopologyEditError::Coordinate(
                        crate::CoordinateValidationError::MissingSourceConformerOrder,
                    ))
                })?;
                conformer
                    .source_remove_atom_position(atom.index(), atom_count)
                    .map_err(|e| E::from(TopologyEditError::Coordinate(e)))?;
                two_d += 1;
            }
            crate::CoordinateDimension::ThreeD => {
                let conformer = coordinates.conformers_3d.get_mut(three_d).ok_or_else(|| {
                    E::from(TopologyEditError::Coordinate(
                        crate::CoordinateValidationError::MissingSourceConformerOrder,
                    ))
                })?;
                conformer
                    .source_remove_atom_position(atom.index(), atom_count)
                    .map_err(|e| E::from(TopologyEditError::Coordinate(e)))?;
                three_d += 1;
            }
        }
    }
    if two_d != coordinates.conformers_2d.len() || three_d != coordinates.conformers_3d.len() {
        return Err(TopologyEditError::Coordinate(
            crate::CoordinateValidationError::MissingSourceConformerOrder,
        )
        .into());
    }
    for (index, bond) in topology.bonds.iter_mut().enumerate() {
        if let Some(value) = bond.prop("_MolFileBondEndPts") {
            let text = text_reader(value)?;
            update_source_bond_endpoints_property(bond, atom.index(), &text).map_err(E::from)?;
        }
        let begin = if bond.begin().index() > atom.index() {
            AtomId::new(bond.begin().index() - 1)
        } else {
            bond.begin()
        };
        let end = if bond.end().index() > atom.index() {
            AtomId::new(bond.end().index() - 1)
        } else {
            bond.end()
        };
        bond.set_endpoints(begin, end);
        bond.set_id_for_construction(BondId::new(index));
        if !bond.stereo_atom_references().is_empty() {
            if bond.stereo_atom_references().contains(&atom) {
                bond.set_stereo_atoms(None);
            } else {
                let refs = bond
                    .stereo_atom_references()
                    .iter()
                    .map(|id| {
                        if id.index() > atom.index() {
                            AtomId::new(id.index() - 1)
                        } else {
                            *id
                        }
                    })
                    .collect();
                bond.set_source_stereo_atom_references(refs);
            }
        }
    }
    crate::sgroup::remove_source_groups_referencing_atom(
        &mut topology.substance_groups,
        atom,
        uint_reader,
    )?;
    for group in &mut topology.stereo_groups {
        group.remove_atom(atom);
        group.shift_source_atom_ids_after_removed(atom);
    }
    topology.stereo_groups.retain(|g| !g.is_empty());
    if clear_props {
        clear_source_computed_properties(topology, properties, true, Some(reset_ring))?;
    }
    topology.atoms.remove(atom.index());
    topology.adjacency = AdjacencyList::from_topology(topology.atoms.len(), &topology.bonds);
    Ok(())
}

#[cfg(test)]
mod source_remove_atom_complete_tests {
    use super::*;
    use crate::{
        BondOrder, Conformer2D, Conformer3D, CoordinateBlock, CoordinateDimension, Element,
        MoleculeProperties, PropertyValue, StereoGroupKind, SubstanceGroupId, SubstanceGroupKind,
    };
    use std::{cell::Cell, collections::BTreeMap};
    fn topology() -> TopologyBlock {
        let atoms = (0..4)
            .map(|i| Atom::from_spec(AtomId::new(i), AtomSpec::new(Element::C)))
            .collect();
        let bonds = (0..3)
            .map(|i| {
                Bond::from_spec(
                    BondId::new(i),
                    BondSpec::new(AtomId::new(i), AtomId::new(i + 1), BondOrder::Single),
                )
            })
            .collect();
        TopologyBlock::try_from_parts(atoms, bonds, vec![], vec![]).unwrap()
    }
    fn uint(value: &PropertyValue) -> Result<u32, TopologyEditError> {
        value
            .as_uint()
            .map_err(|_| TopologyEditError::SourceGroupParentCycle {
                group: SubstanceGroupId::new(999),
            })
    }
    fn text(value: &PropertyValue) -> Result<crate::PropertyText, TopologyEditError> {
        Ok(value
            .as_string()
            .expect("fixture source reader uses exact StringTag")
            .clone())
    }
    fn coords() -> CoordinateBlock {
        CoordinateBlock {
            conformers_3d: vec![
                Conformer3D::new(
                    27,
                    vec![
                        [-0.0, 1.0, 2.0],
                        [10.0, 11.0, 12.0],
                        [20.0, 21.0, 22.0],
                        [30.0, 31.0, 32.0],
                    ],
                    true,
                )
                .with_prop("retained", "value"),
            ],
            ..Default::default()
        }
    }
    #[test]
    fn source_remove_atom_batch_schedules_incident_bonds_before_atom_without_other_effects() {
        let mut topology = topology();
        let before = topology.clone();
        let mut atoms = vec![];
        let mut bonds = vec![];
        remove_atom_source::<TopologyEditError>(
            &mut topology,
            AtomId::new(1),
            false,
            SourceAtomRemovalState::Batch {
                atoms: &mut atoms,
                bonds: &mut bonds,
            },
        )
        .unwrap();
        assert_eq!(atoms, [false, true, false, false]);
        assert_eq!(bonds, [true, true, false]);
        assert_eq!(topology, before);
        let mut edit = topology.into_batch_edit().unwrap();
        edit.remove_atom(AtomId::new(1)).unwrap();
        assert_eq!(edit.remove_bonds, [true, true, false]);
        assert!(edit.remove_atoms[1]);
    }
    #[test]
    fn source_remove_atom_immediate_updates_all_references_coordinates_endpts_and_clear_flag() {
        for clear in [false, true] {
            let mut topology = topology();
            topology.atoms[3]
                .set_computed_prop("atom-computed", PropertyValue::UInt(5))
                .unwrap();
            topology.bonds[2]
                .set_prop(
                    "_MolFileBondEndPts",
                    PropertyValue::String("(4 1 2 4 2)".into()),
                )
                .unwrap();
            topology.bonds[2]
                .set_computed_prop("bond-computed", PropertyValue::UInt(6))
                .unwrap();
            topology.substance_groups.push(
                SubstanceGroup::new(SubstanceGroupId::new(0), SubstanceGroupKind::Data)
                    .with_atoms(vec![AtomId::new(3), AtomId::new(2)]),
            );
            topology.stereo_groups.push(
                StereoGroup::new(
                    StereoGroupKind::And,
                    vec![AtomId::new(1), AtomId::new(3)],
                    vec![BondId::new(2)],
                )
                .expect("valid distinct stereo members"),
            );
            let mut coordinates = coords();
            let conf_before = coordinates.conformers_3d[0].clone();
            let mut properties = MoleculeProperties::default();
            properties
                .set_computed_prop("mol-computed", PropertyValue::UInt(7))
                .unwrap();
            let mut atoms =
                BTreeMap::from([(10, vec![AtomId::new(1)]), (11, vec![AtomId::new(3)])]);
            let mut bonds = BTreeMap::from([(12, vec![BondId::new(2)])]);
            let resets = Cell::new(0);
            remove_atom_source(
                &mut topology,
                AtomId::new(1),
                clear,
                SourceAtomRemovalState::Immediate {
                    atom_bookmarks: &mut atoms,
                    bond_bookmarks: &mut bonds,
                    coordinates: &mut coordinates,
                    properties: &mut properties,
                    reset_ring: &mut || resets.set(resets.get() + 1),
                    uint_reader: &mut uint,
                    text_reader: &mut text,
                },
            )
            .unwrap();
            assert_eq!(topology.atoms.len(), 3);
            assert_eq!(topology.atoms[2].id(), AtomId::new(2));
            assert_eq!(
                (topology.bonds[0].begin(), topology.bonds[0].end()),
                (AtomId::new(1), AtomId::new(2))
            );
            assert_eq!(
                topology.bonds[0].prop("_MolFileBondEndPts"),
                Some(&PropertyValue::String("(3 1 3 2)".into()))
            );
            assert_eq!(
                topology.substance_groups[0].atoms(),
                [AtomId::new(2), AtomId::new(1)]
            );
            assert_eq!(topology.stereo_groups[0].atoms(), [AtomId::new(2)]);
            assert_eq!(topology.stereo_groups[0].bonds(), [BondId::new(0)]);
            assert_eq!(atoms, BTreeMap::from([(11, vec![AtomId::new(2)])]));
            assert_eq!(bonds, BTreeMap::from([(12, vec![BondId::new(0)])]));
            let conf = &coordinates.conformers_3d[0];
            assert_eq!(conf.id(), 27);
            assert_eq!(conf.props(), conf_before.props());
            assert_eq!(
                conf.coordinates(),
                [
                    conf_before.coordinates()[0],
                    conf_before.coordinates()[2],
                    conf_before.coordinates()[3]
                ]
            );
            assert_eq!(conf.coordinates()[0][0].to_bits(), (-0.0f64).to_bits());
            assert_eq!(resets.get(), if clear { 3 } else { 2 });
            assert_eq!(topology.atoms[2].prop("atom-computed").is_some(), !clear);
            assert_eq!(topology.bonds[0].prop("bond-computed").is_some(), !clear);
            assert_eq!(
                properties.props().contains_key(b"mol-computed".as_slice()),
                !clear
            );
            topology.validate().unwrap();
        }
    }
    #[test]
    fn source_remove_atom_endpts_failure_preserves_prior_bond_bookmark_index_and_coordinate_effects()
     {
        let mut topology = topology();
        topology.bonds[2]
            .set_prop(
                "_MolFileBondEndPts",
                PropertyValue::String("(bad 4)".into()),
            )
            .unwrap();
        let mut coordinates = coords();
        let mut properties = MoleculeProperties::default();
        let mut atoms = BTreeMap::from([(1, vec![AtomId::new(1)]), (2, vec![AtomId::new(3)])]);
        let mut bonds = BTreeMap::new();
        let resets = Cell::new(0);
        let error = remove_atom_source(
            &mut topology,
            AtomId::new(1),
            true,
            SourceAtomRemovalState::Immediate {
                atom_bookmarks: &mut atoms,
                bond_bookmarks: &mut bonds,
                coordinates: &mut coordinates,
                properties: &mut properties,
                reset_ring: &mut || resets.set(resets.get() + 1),
                uint_reader: &mut uint,
                text_reader: &mut text,
            },
        )
        .unwrap_err();
        assert!(matches!(
            error,
            TopologyEditError::BondEndPointsParse { token_index: 0, .. }
        ));
        assert_eq!(topology.atoms.len(), 4);
        assert_eq!(topology.atoms[3].id(), AtomId::new(2));
        assert_eq!(coordinates.conformers_3d[0].coordinates().len(), 3);
        assert_eq!(topology.bonds.len(), 1);
        assert_eq!(
            (topology.bonds[0].begin(), topology.bonds[0].end()),
            (AtomId::new(2), AtomId::new(3))
        );
        assert_eq!(atoms, BTreeMap::from([(2, vec![AtomId::new(2)])]));
        assert_eq!(resets.get(), 2);
    }
    #[test]
    fn source_remove_atom_mixed_conformer_failure_keeps_prior_source_order_erase() {
        let mut topology = topology();
        let mut coordinates = coords();
        coordinates.conformers_2d.push(Conformer2D::new(8, vec![]));
        coordinates.source_conformer_order =
            Some(vec![CoordinateDimension::ThreeD, CoordinateDimension::TwoD]);
        let mut properties = MoleculeProperties::default();
        let mut atoms = BTreeMap::new();
        let mut bonds = BTreeMap::new();
        let error = remove_atom_source(
            &mut topology,
            AtomId::new(0),
            false,
            SourceAtomRemovalState::Immediate {
                atom_bookmarks: &mut atoms,
                bond_bookmarks: &mut bonds,
                coordinates: &mut coordinates,
                properties: &mut properties,
                reset_ring: &mut || {},
                uint_reader: &mut uint,
                text_reader: &mut text,
            },
        )
        .unwrap_err();
        assert!(matches!(
            error,
            TopologyEditError::Coordinate(crate::CoordinateValidationError::RowCount {
                dimension: "2D",
                ..
            })
        ));
        assert_eq!(coordinates.conformers_3d[0].coordinates().len(), 3);
        assert_eq!(coordinates.conformers_2d[0].coordinates().len(), 0);
        assert_eq!(topology.atoms.len(), 4);
    }
    #[test]
    fn source_remove_atom_computed_property_error_follows_reference_edits_before_vertex_erasure() {
        let mut topology = topology();
        let mut coordinates = coords();
        let mut properties = MoleculeProperties::default();
        properties
            .set_prop("__computedProps", PropertyValue::UInt(99))
            .unwrap();
        let mut atoms = BTreeMap::new();
        let mut bonds = BTreeMap::new();
        let resets = Cell::new(0);
        let error = remove_atom_source_index(
            &mut topology,
            AtomId::new(1),
            SourceAtomRemovalState::Immediate {
                atom_bookmarks: &mut atoms,
                bond_bookmarks: &mut bonds,
                coordinates: &mut coordinates,
                properties: &mut properties,
                reset_ring: &mut || resets.set(resets.get() + 1),
                uint_reader: &mut uint,
                text_reader: &mut text,
            },
        )
        .unwrap_err();
        assert!(matches!(error, TopologyEditError::MoleculeProperty(_)));
        assert_eq!(resets.get(), 3);
        assert_eq!(topology.atoms.len(), 4);
        assert_eq!(
            (topology.bonds[0].begin(), topology.bonds[0].end()),
            (AtomId::new(1), AtomId::new(2))
        );
        assert_eq!(coordinates.conformers_3d[0].coordinates().len(), 3);
    }
}

#[cfg(test)]
mod source_bookmark_erasure_tests {
    use super::*;
    use crate::{BondOrder, CoordinateBlock, Element, MoleculeProperties, PropertyValue};
    use std::collections::BTreeMap;
    fn topology() -> TopologyBlock {
        TopologyBlock::try_from_parts(
            vec![
                Atom::from_spec(AtomId::new(0), AtomSpec::new(Element::C)),
                Atom::from_spec(AtomId::new(1), AtomSpec::new(Element::N)),
            ],
            vec![Bond::from_spec(
                BondId::new(0),
                BondSpec::new(AtomId::new(0), AtomId::new(1), BondOrder::Single),
            )],
            vec![],
            vec![],
        )
        .unwrap()
    }
    fn uint(value: &PropertyValue) -> Result<u32, TopologyEditError> {
        Ok(value.as_uint().unwrap())
    }
    #[test]
    fn source_bookmark_erasure_retains_unrelated_empty_map_entries() {
        // ROMol.cpp clear{Atom,Bond}Bookmark is reached only when the outer
        // remove loop finds the target pointer. An unrelated empty list stays.
        let mut topology = topology();
        let mut bonds = BTreeMap::from([(1, vec![]), (2, vec![BondId::new(0)])]);
        remove_bond_source(
            &mut topology,
            AtomId::new(0),
            AtomId::new(1),
            SourceBondRemovalState::Immediate {
                bookmarks: &mut bonds,
                reset_ring: &mut || {},
                uint_reader: &mut uint,
            },
        )
        .unwrap();
        assert_eq!(bonds, BTreeMap::from([(1, vec![])]));
        let mut atoms = BTreeMap::from([(3, vec![]), (4, vec![AtomId::new(0)])]);
        let mut coordinates = CoordinateBlock::default();
        let mut props = MoleculeProperties::default();
        remove_atom_source(
            &mut topology,
            AtomId::new(0),
            false,
            SourceAtomRemovalState::Immediate {
                atom_bookmarks: &mut atoms,
                bond_bookmarks: &mut bonds,
                coordinates: &mut coordinates,
                properties: &mut props,
                reset_ring: &mut || {},
                uint_reader: &mut uint,
                text_reader: &mut |v| Ok(v.as_string().unwrap().clone()),
            },
        )
        .unwrap();
        assert_eq!(atoms, BTreeMap::from([(3, vec![])]));
        assert_eq!(bonds, BTreeMap::from([(1, vec![])]));
    }
}

/// Source graph adjacency borrows original CSR plus ordered appended edges.
/// Both views are detached canonical rows, without live/cache authority.
#[doc(hidden)]
pub struct SourceBondNeighbors<'a> {
    pub original: Option<&'a AdjacencyList>,
    pub appended: &'a mut [Vec<crate::NeighborRef>],
}
#[doc(hidden)]
pub struct SourceBondBatchMasks<'a> {
    pub atoms: &'a [bool],
    pub bonds: &'a mut Vec<bool>,
}

/// Canonical source order overload reused by detached construction consumers.
#[doc(hidden)]
pub fn add_source_bond_order(
    atoms: &mut [Atom],
    bonds: &mut Vec<Bond>,
    neighbors: SourceBondNeighbors<'_>,
    pending: Option<SourceBondBatchMasks<'_>>,
    begin: AtomId,
    end: AtomId,
    order: crate::BondOrder,
) -> Result<usize, TopologyEditError> {
    // RDKit❗✔️: unsigned int RWMol::addBond(unsigned int atomIdx1, unsigned int atomIdx2,
    // RDKit❗✔️:                             Bond::BondType bondType) {
    // RDKit❗✔️:   // if the atom indices are bad, the next two calls will catch that.
    // RDKit❗✔️:   auto beginAtom = getAtomWithIdx(atomIdx1);
    // RDKit❗✔️:   auto endAtom = getAtomWithIdx(atomIdx2);
    // RDKit❗✔️:   PRECONDITION(atomIdx1 != atomIdx2, "attempt to add self-bond");
    // RDKit❗✔️:   PRECONDITION(!(boost::edge(atomIdx1, atomIdx2, d_graph).second),
    // RDKit❗✔️:                "bond already exists");
    // RDKit❗✔️:
    // RDKit❗✔️:   auto *b = new Bond(bondType);
    // RDKit❗✔️:   b->setOwningMol(this);
    // RDKit❗✔️:   if (bondType == Bond::AROMATIC) {
    // RDKit❗✔️:     b->setIsAromatic(1);
    // RDKit❗✔️:     //
    // RDKit❗✔️:     // assume that aromatic bonds connect aromatic atoms
    // RDKit❗✔️:     //   This is relevant for file formats like MOL, where there
    // RDKit❗✔️:     //   is no such thing as an aromatic atom, but bonds can be
    // RDKit❗✔️:     //   marked aromatic.
    // RDKit❗✔️:     //
    // RDKit❗✔️:     beginAtom->setIsAromatic(1);
    // RDKit❗✔️:     endAtom->setIsAromatic(1);
    // RDKit❗✔️:   }
    // RDKit❗✔️:   auto [which, ok] = boost::add_edge(atomIdx1, atomIdx2, d_graph);
    // RDKit❗✔️:   d_graph[which] = b;
    // RDKit❗✔️:   ++numBonds;
    // RDKit❗✔️:   b->setIdx(numBonds - 1);
    // RDKit❗✔️:   b->setBeginAtomIdx(atomIdx1);
    // RDKit❗✔️:   b->setEndAtomIdx(atomIdx2);
    // RDKit❗✔️:
    // RDKit❗✔️:   // the valence values on the begin and end atoms need to be updated:
    // RDKit❗✔️:   beginAtom->clearPropertyCache();
    // RDKit❗✔️:   endAtom->clearPropertyCache();
    // RDKit❗✔️:
    // RDKit❗✔️:   // we're in a batch edit, and at least one of the bond ends is scheduled
    // RDKit❗✔️:   // for deletion, so mark the new bond for deletion too:
    // RDKit❗✔️:   if (dp_delAtoms &&
    // RDKit❗✔️:       ((atomIdx1 < dp_delAtoms->size() && dp_delAtoms->test(atomIdx1)) ||
    // RDKit❗✔️:        (atomIdx2 < dp_delAtoms->size() && dp_delAtoms->test(atomIdx2)))) {
    // RDKit❗✔️:     if (dp_delBonds->size() < numBonds) {
    // RDKit❗✔️:       dp_delBonds->resize(numBonds);
    // RDKit❗✔️:     }
    // RDKit❗✔️:     dp_delBonds->set(numBonds - 1);
    // RDKit❗✔️:   }
    // RDKit❗✔️:
    // RDKit❗✔️:   return numBonds;
    // RDKit❗✔️: }
    // Native pointers/unsigned32 graph width/allocator failures remain gaps.
    // Cost: O(degree) borrowed lookup and amortized row/two-neighbor appends;
    // two scalar source cache resets and conditional mask resize, no graph clone.
    for atom in [begin, end] {
        if atom.index() >= atoms.len() {
            return Err(TopologyEditError::AtomOutOfRange {
                atom,
                atom_count: atoms.len(),
            });
        }
    }
    if begin == end {
        return Err(TopologyEditError::InvalidResult(
            TopologyValidationError::SelfLoopBond {
                bond: BondId::new(bonds.len()),
                atom: begin,
            },
        ));
    }
    if neighbors.appended.len() != atoms.len() {
        return Err(TopologyEditError::InvalidSource(
            TopologyValidationError::AdjacencyMismatch,
        ));
    }
    let duplicate = neighbors.original.is_some_and(|adjacency| {
        adjacency
            .neighbors_of(begin.index())
            .iter()
            .any(|n| n.atom_index == end.index())
    }) || neighbors.appended[begin.index()]
        .iter()
        .any(|n| n.atom_index == end.index());
    if duplicate {
        return Err(TopologyEditError::DuplicateBond { begin, end });
    }
    let aromatic = order == crate::BondOrder::Aromatic;
    let id = BondId::new(bonds.len());
    let bond = Bond::from_spec(id, BondSpec::new(begin, end, order).with_aromatic(aromatic));
    if aromatic {
        atoms[begin.index()].set_aromatic(true);
        atoms[end.index()].set_aromatic(true);
    }
    bonds.push(bond);
    neighbors.appended[begin.index()].push(crate::NeighborRef {
        atom_index: end.index(),
        bond: id,
    });
    neighbors.appended[end.index()].push(crate::NeighborRef {
        atom_index: begin.index(),
        bond: id,
    });
    // RDKit✔️✔️: void Atom::clearPropertyCache() {
    // RDKit✔️✔️:   d_explicitValence = -1;
    // RDKit✔️✔️:   d_implicitValence = -1;
    // RDKit✔️✔️: }
    atoms[begin.index()].set_source_valence_facts(crate::SourceAtomValenceFacts::UNINITIALIZED);
    atoms[end.index()].set_source_valence_facts(crate::SourceAtomValenceFacts::UNINITIALIZED);
    if let Some(masks) = pending {
        if masks.atoms.get(begin.index()).copied().unwrap_or(false)
            || masks.atoms.get(end.index()).copied().unwrap_or(false)
        {
            if masks.bonds.len() < bonds.len() {
                masks.bonds.resize(bonds.len(), false);
            }
            masks.bonds[bonds.len() - 1] = true;
        }
    }
    Ok(bonds.len())
}

#[cfg(test)]
mod source_add_bond_order_complete_tests {
    use super::*;
    use crate::{BondOrder, Element, SourceAtomValenceFacts};
    fn atoms() -> Vec<Atom> {
        (0..3)
            .map(|i| {
                let mut atom = Atom::from_spec(AtomId::new(i), AtomSpec::new(Element::C));
                atom.set_source_valence_facts(SourceAtomValenceFacts {
                    explicit_valence: 3,
                    implicit_valence: 1,
                });
                atom
            })
            .collect()
    }
    #[test]
    fn source_add_bond_order_aromatic_endpoint_flags_cache_reset_and_count() {
        for order in [
            BondOrder::Single,
            BondOrder::Aromatic,
            BondOrder::Unspecified,
            BondOrder::Other,
            BondOrder::DativeLeft,
        ] {
            let mut atoms = atoms();
            let mut bonds = Vec::new();
            let mut neighbors = vec![vec![]; 3];
            let count = add_source_bond_order(
                &mut atoms,
                &mut bonds,
                SourceBondNeighbors {
                    original: None,
                    appended: &mut neighbors,
                },
                None,
                AtomId::new(1),
                AtomId::new(2),
                order,
            )
            .unwrap();
            assert_eq!(count, 1);
            assert_eq!(bonds[0].id(), BondId::new(0));
            assert_eq!(bonds[0].order(), order);
            assert_eq!(bonds[0].is_aromatic(), order == BondOrder::Aromatic);
            assert_eq!(atoms[1].is_aromatic(), order == BondOrder::Aromatic);
            assert_eq!(atoms[2].is_aromatic(), order == BondOrder::Aromatic);
            assert_eq!(
                atoms[1].source_valence_facts(),
                SourceAtomValenceFacts::UNINITIALIZED
            );
            assert_eq!(
                atoms[2].source_valence_facts(),
                SourceAtomValenceFacts::UNINITIALIZED
            );
            assert_eq!(atoms[0].source_valence_facts().explicit_valence, 3);
            assert_eq!(neighbors[1][0].atom_index, 2);
            assert_eq!(neighbors[2][0].atom_index, 1);
        }
    }
    #[test]
    fn source_add_bond_order_mask_resizes_only_for_a_scheduled_endpoint() {
        let mut atoms = atoms();
        let mut bonds = Vec::new();
        let mut neighbors = vec![vec![]; 3];
        let mut pending = vec![];
        add_source_bond_order(
            &mut atoms,
            &mut bonds,
            SourceBondNeighbors {
                original: None,
                appended: &mut neighbors,
            },
            Some(SourceBondBatchMasks {
                atoms: &[],
                bonds: &mut pending,
            }),
            AtomId::new(0),
            AtomId::new(1),
            BondOrder::Single,
        )
        .unwrap();
        assert!(pending.is_empty()); // native source does not resize unmarked mask
        let count = add_source_bond_order(
            &mut atoms,
            &mut bonds,
            SourceBondNeighbors {
                original: None,
                appended: &mut neighbors,
            },
            Some(SourceBondBatchMasks {
                atoms: &[false, true],
                bonds: &mut pending,
            }),
            AtomId::new(1),
            AtomId::new(2),
            BondOrder::Aromatic,
        )
        .unwrap();
        assert_eq!(count, 2);
        assert_eq!(pending, [false, true]);
        let mut edit = TopologyBlock::try_from_parts(atoms, bonds, vec![], vec![])
            .unwrap()
            .into_batch_edit()
            .unwrap();
        edit.remove_atom(AtomId::new(2)).unwrap();
        assert_eq!(
            edit.add_bond_by_order(AtomId::new(0), AtomId::new(2), BondOrder::Double)
                .unwrap(),
            3
        );
        assert!(edit.remove_bonds[2]);
        assert_eq!(
            edit.bond_between_atoms(AtomId::new(0), AtomId::new(2))
                .unwrap(),
            Some(BondId::new(2))
        );
    }
    #[test]
    fn source_add_bond_order_endpoint_self_duplicate_checks_precede_any_mutation() {
        let mut atoms = atoms();
        let mut bonds = Vec::new();
        let mut neighbors = vec![vec![]; 3];
        let error = add_source_bond_order(
            &mut atoms,
            &mut bonds,
            SourceBondNeighbors {
                original: None,
                appended: &mut neighbors,
            },
            None,
            AtomId::new(9),
            AtomId::new(8),
            BondOrder::Aromatic,
        )
        .unwrap_err();
        assert!(
            matches!(error, TopologyEditError::AtomOutOfRange { atom, .. } if atom == AtomId::new(9))
        );
        let error = add_source_bond_order(
            &mut atoms,
            &mut bonds,
            SourceBondNeighbors {
                original: None,
                appended: &mut neighbors,
            },
            None,
            AtomId::new(1),
            AtomId::new(1),
            BondOrder::Aromatic,
        )
        .unwrap_err();
        assert!(matches!(
            error,
            TopologyEditError::InvalidResult(TopologyValidationError::SelfLoopBond { .. })
        ));
        assert!(bonds.is_empty());
        assert!(
            atoms
                .iter()
                .all(|a| !a.is_aromatic() && a.source_valence_facts().explicit_valence == 3)
        );
        add_source_bond_order(
            &mut atoms,
            &mut bonds,
            SourceBondNeighbors {
                original: None,
                appended: &mut neighbors,
            },
            None,
            AtomId::new(0),
            AtomId::new(1),
            BondOrder::Single,
        )
        .unwrap();
        let before_atoms = atoms.clone();
        let before_bonds = bonds.clone();
        let before_neighbors = neighbors.clone();
        let error = add_source_bond_order(
            &mut atoms,
            &mut bonds,
            SourceBondNeighbors {
                original: None,
                appended: &mut neighbors,
            },
            None,
            AtomId::new(1),
            AtomId::new(0),
            BondOrder::Aromatic,
        )
        .unwrap_err();
        assert!(matches!(error, TopologyEditError::DuplicateBond { .. }));
        assert_eq!(atoms, before_atoms);
        assert_eq!(bonds, before_bonds);
        assert_eq!(neighbors, before_neighbors);
    }
}

#[cfg(test)]
mod source_add_bond_atom_overload_tests {
    use super::*;
    use crate::{BondOrder, Element, SourceAtomValenceFacts};
    #[test]
    fn source_atom_bond_overload_null_precedes_bad_id_and_foreign_atom_identity_is_only_its_index()
    {
        let mut edit = TopologyBlock::try_from_parts(
            vec![
                Atom::from_spec(AtomId::new(0), AtomSpec::new(Element::C)),
                Atom::from_spec(AtomId::new(1), AtomSpec::new(Element::C)),
            ],
            vec![],
            vec![],
            vec![],
        )
        .unwrap()
        .into_batch_edit()
        .unwrap();
        let bad = Atom::from_spec(AtomId::new(9), AtomSpec::new(Element::N));
        for (begin, end) in [(None, None), (Some(&bad), None), (None, Some(&bad))] {
            assert!(matches!(
                edit.add_bond_from_atoms(begin, end, BondOrder::Aromatic),
                Err(TopologyEditError::NullBondAtom)
            ));
            assert!(edit.working.bonds.is_empty());
            assert!(!edit.working.atoms[0].is_aromatic());
        }
        let mut foreign_begin = Atom::from_spec(AtomId::new(0), AtomSpec::new(Element::N));
        let mut foreign_end = Atom::from_spec(AtomId::new(1), AtomSpec::new(Element::O));
        foreign_begin.set_source_valence_facts(SourceAtomValenceFacts {
            explicit_valence: 3,
            implicit_valence: 1,
        });
        foreign_end.set_source_valence_facts(SourceAtomValenceFacts {
            explicit_valence: 2,
            implicit_valence: 0,
        });
        let original_begin = foreign_begin.clone();
        let original_end = foreign_end.clone();
        assert_eq!(
            edit.add_bond_from_atoms(
                Some(&foreign_begin),
                Some(&foreign_end),
                BondOrder::Aromatic
            )
            .unwrap(),
            1
        );
        assert_eq!(foreign_begin, original_begin);
        assert_eq!(foreign_end, original_end);
        assert!(edit.working.atoms[0].is_aromatic() && edit.working.atoms[1].is_aromatic());
        assert_eq!(edit.working.atoms[0].element(), Element::C);
        assert_eq!(edit.working.atoms[1].element(), Element::C);
        assert_eq!(edit.working.bonds[0].id(), BondId::new(0));
    }
}

fn clear_source_bond_bookmarks(
    bookmarks: &mut std::collections::BTreeMap<i32, Vec<BondId>>,
    bond_id: BondId,
) {
    // RDKit❗✔️: if (std::find(bonds.begin(), bonds.end(), bnd) != bonds.end()) {
    // RDKit❗✔️:   clearBondBookmark(tmpI->first, bnd);
    // RDKit❗✔️: }
    bookmarks.retain(|_, values| {
        if let Some(position) = values.iter().position(|b| *b == bond_id) {
            values.remove(position);
            !values.is_empty()
        } else {
            true
        }
    });
}

fn clear_source_incident_bond_stereo(
    topology: &mut TopologyBlock,
    index: usize,
) -> Result<(), TopologyEditError> {
    // RDKit❗❌:     auto beginAtm = bnd->getBeginAtom();
    // RDKit❗❌:     auto endAtm = bnd->getEndAtom();
    // RDKit❗❌:     std::vector<std::vector<Atom *>> bond_atoms = {{beginAtm, endAtm},
    // RDKit❗❌:                                                    {endAtm, beginAtm}};
    // RDKit❗❌:     for (const auto &atoms : bond_atoms) {
    // RDKit❗❌:       for (auto obnd : atomBonds(atoms[0])) {
    // RDKit❗❌:         if (obnd == bnd) {
    // RDKit❗❌:           continue;
    // RDKit❗❌:         }
    // RDKit❗❌:         if (std::find(obnd->getStereoAtoms().begin(),
    // RDKit❗❌:                       obnd->getStereoAtoms().end(),
    // RDKit❗❌:                       atoms[1]->getIdx()) != obnd->getStereoAtoms().end()) {
    // RDKit❗❌:           // github #6900 if we remove stereo atoms we need to remove
    // RDKit❗❌:           //  the CIS and or TRANS since this requires stereo atoms
    // RDKit❗❌:           if (obnd->getStereo() == Bond::BondStereo::STEREOCIS ||
    // RDKit❗❌:               obnd->getStereo() == Bond::BondStereo::STEREOTRANS) {
    // RDKit❗❌:             obnd->setStereo(Bond::BondStereo::STEREONONE);
    // RDKit❗❌:           }
    // RDKit❗❌:           obnd->getStereoAtoms().clear();
    // RDKit❗❌:         }
    // RDKit❗❌:       }
    // RDKit❗❌:     }
    // RDKit❗❌:
    let source_begin = topology.bonds[index].begin();
    let source_end = topology.bonds[index].end();
    for (endpoint, removed_neighbor) in [(source_begin, source_end), (source_end, source_begin)] {
        for (other_index, other) in topology.bonds.iter_mut().enumerate() {
            if other_index == index || (other.begin() != endpoint && other.end() != endpoint) {
                continue;
            }
            if other.stereo_atom_references().contains(&removed_neighbor) {
                if matches!(other.stereo(), BondStereo::Cis | BondStereo::Trans) {
                    other
                        .set_stereo(BondStereo::None)
                        .map_err(TopologyEditError::InvalidBond)?;
                }
                other.set_stereo_atoms(None);
            }
        }
    }
    Ok(())
}

/// Source bond-only deletion reused by detached construction owners.
#[doc(hidden)]
pub fn batch_remove_bonds_source<E: From<TopologyEditError>>(
    topology: &mut TopologyBlock,
    pending: Option<&[bool]>,
    mut bookmarks: Option<&mut std::collections::BTreeMap<i32, Vec<BondId>>>,
    uint_reader: &mut (impl FnMut(&crate::PropertyValue) -> Result<u32, E> + ?Sized),
) -> Result<(), E> {
    // RDKit❗❌: void RWMol::batchRemoveBonds() {
    // RDKit❗❌:   if (!dp_delBonds || dp_delBonds->none()) {
    // RDKit❗❌:     return;
    // RDKit❗❌:   }
    // RDKit❗❌:
    // RDKit❗❌:   auto &delBonds = *dp_delBonds;
    // RDKit❗❌:   unsigned int min_idx = getNumBonds();
    // RDKit❗❌:   for (unsigned int i = rdcast<unsigned int>(delBonds.size()); i > 0; --i) {
    // RDKit❗❌:     if (!delBonds[i - 1]) {
    // RDKit❗❌:       continue;
    // RDKit❗❌:     }
    // RDKit❗❌:     unsigned int idx = rdcast<unsigned int>(i - 1);
    // RDKit❗❌:     Bond *bnd = getBondWithIdx(idx);
    // RDKit❗❌:     if (!bnd) {
    // RDKit❗❌:       continue;
    // RDKit❗❌:     }
    // RDKit❗❌:
    // RDKit❗❌:     min_idx = idx;
    // RDKit❗❌:     // remove any bookmarks which point to this bond:
    // RDKit❗❌:     BOND_BOOKMARK_MAP *marks = getBondBookmarks();
    // RDKit❗❌:     auto markI = marks->begin();
    // RDKit❗❌:     while (markI != marks->end()) {
    // RDKit❗❌:       BOND_PTR_LIST &bonds = markI->second;
    // RDKit❗❌:       // we need to copy the iterator then increment it, because the
    // RDKit❗❌:       // deletion we're going to do in clearBondBookmark will invalidate
    // RDKit❗❌:       // it.
    // RDKit❗❌:       auto tmpI = markI;
    // RDKit❗❌:       ++markI;
    // RDKit❗❌:       if (std::find(bonds.begin(), bonds.end(), bnd) != bonds.end()) {
    // RDKit❗❌:         clearBondBookmark(tmpI->first, bnd);
    // RDKit❗❌:       }
    // RDKit❗❌:     }
    // RDKit❗❌:
    // RDKit❗❌:     // loop over neighboring double bonds and remove their stereo atom
    // RDKit❗❌:     //  information. This is definitely now invalid (was github issue 8)
    // RDKit❗❌:     auto beginAtm = bnd->getBeginAtom();
    // RDKit❗❌:     auto endAtm = bnd->getEndAtom();
    // RDKit❗❌:     std::vector<std::vector<Atom *>> bond_atoms = {{beginAtm, endAtm},
    // RDKit❗❌:                                                    {endAtm, beginAtm}};
    // RDKit❗❌:     for (const auto &atoms : bond_atoms) {
    // RDKit❗❌:       for (auto obnd : atomBonds(atoms[0])) {
    // RDKit❗❌:         if (obnd == bnd) {
    // RDKit❗❌:           continue;
    // RDKit❗❌:         }
    // RDKit❗❌:         if (std::find(obnd->getStereoAtoms().begin(),
    // RDKit❗❌:                       obnd->getStereoAtoms().end(),
    // RDKit❗❌:                       atoms[1]->getIdx()) != obnd->getStereoAtoms().end()) {
    // RDKit❗❌:           // github #6900 if we remove stereo atoms we need to remove
    // RDKit❗❌:           //  the CIS and or TRANS since this requires stereo atoms
    // RDKit❗❌:           if (obnd->getStereo() == Bond::BondStereo::STEREOCIS ||
    // RDKit❗❌:               obnd->getStereo() == Bond::BondStereo::STEREOTRANS) {
    // RDKit❗❌:             obnd->setStereo(Bond::BondStereo::STEREONONE);
    // RDKit❗❌:           }
    // RDKit❗❌:           obnd->getStereoAtoms().clear();
    // RDKit❗❌:         }
    // RDKit❗❌:       }
    // RDKit❗❌:     }
    // RDKit❗❌:
    // RDKit❗❌:     removeSubstanceGroupsReferencingBond(*this, idx);
    // RDKit❗❌:     // Remove this bond from any stereo group
    // RDKit❗❌:     removeBondFromGroups(bnd, d_stereo_groups);
    // RDKit❗❌:
    // RDKit❗❌:     bnd->setOwningMol(nullptr);
    // RDKit❗❌:
    // RDKit❗❌:     auto vd1 = boost::vertex(bnd->getBeginAtomIdx(), d_graph);
    // RDKit❗❌:     auto vd2 = boost::vertex(bnd->getEndAtomIdx(), d_graph);
    // RDKit❗❌:     boost::remove_edge(vd1, vd2, d_graph);
    // RDKit❗❌:     delete bnd;
    // RDKit❗❌:     --numBonds;
    // RDKit❗❌:   }
    // RDKit❗❌:
    // RDKit❗❌:   // loop over all bonds with higher indices than the minimum modified and
    // RDKit❗❌:   // update their indices
    // RDKit❗❌:   auto [firstB, lastB] = this->getEdges();
    // RDKit❗❌:   unsigned int next_idx = min_idx;
    // RDKit❗❌:   while (firstB != lastB) {
    // RDKit❗❌:     Bond *bond = (*this)[*firstB];
    // RDKit❗❌:     if (bond->getIdx() > min_idx) {
    // RDKit❗❌:       bond->setIdx(next_idx++);
    // RDKit❗❌:     }
    // RDKit❗❌:     ++firstB;
    // RDKit❗❌:   }
    // RDKit❗❌: }
    // Native pointer/null/width/allocator states are unmodeled. This actual
    // source batch method is separate from the local structural mapping finish;
    // native consumer integration remains at the whole source gate.
    // Cost: per-delete canonical row stereo scans/Vec erase O(E), plus one dense
    // ID alias map/CSR rebuild. These add source-absent work; no graph clone.
    let Some(pending) = pending else {
        return Ok(());
    };
    if !pending.iter().any(|removed| *removed) {
        return Ok(());
    }
    let original_count = topology.bonds.len();
    let mut minimum = original_count;
    for index in (0..pending.len()).rev() {
        if !pending[index] {
            continue;
        }
        let bond = topology.bonds.get(index).ok_or_else(|| {
            E::from(TopologyEditError::BondOutOfRange {
                bond: BondId::new(index),
                bond_count: topology.bonds.len(),
            })
        })?;
        let old_id = bond.id();
        minimum = index;
        if let Some(bookmarks) = bookmarks.as_deref_mut() {
            clear_source_bond_bookmarks(bookmarks, old_id);
        }
        clear_source_incident_bond_stereo(topology, index).map_err(E::from)?;
        crate::sgroup::remove_source_groups_referencing_bond(
            &mut topology.substance_groups,
            BondId::new(index),
            uint_reader,
        )?;
        for group in &mut topology.stereo_groups {
            group.remove_bond(old_id);
        }
        topology.stereo_groups.retain(|g| !g.is_empty());
        topology.bonds.remove(index);
    }
    let mut aliases = vec![None; original_count];
    let mut next = minimum;
    for bond in &mut topology.bonds {
        let old = bond.id();
        if old.index() > minimum {
            bond.set_id_for_construction(BondId::new(next));
            next += 1;
        }
        let Some(alias) = aliases.get_mut(old.index()) else {
            return Err(TopologyEditError::InvalidSource(
                TopologyValidationError::BondIdMismatch {
                    position: old.index(),
                    id: old,
                },
            )
            .into());
        };
        *alias = Some(bond.id());
    }
    // Native pointer aliases observe final setIdx automatically. Detached IDs
    // track retained objects. Surviving references to deleted pointers are a
    // native dangling-state gap and are not silently removed or reinterpreted.
    if let Some(bookmarks) = bookmarks {
        for values in bookmarks.values_mut() {
            for id in values {
                if let Some(Some(new_id)) = aliases.get(id.index()) {
                    *id = *new_id;
                }
            }
        }
    }
    for group in &mut topology.stereo_groups {
        group.remap_source_retained_bond_ids(&aliases);
    }
    topology.adjacency = AdjacencyList::from_topology(topology.atoms.len(), &topology.bonds);
    Ok(())
}

#[cfg(test)]
mod source_batch_remove_bonds_complete_tests {
    use super::*;
    use crate::{
        BondOrder, Element, PropertyValue, StereoGroupKind, SubstanceGroupId, SubstanceGroupKind,
    };
    use std::collections::BTreeMap;
    fn topology() -> TopologyBlock {
        let atoms = (0..4)
            .map(|i| Atom::from_spec(AtomId::new(i), AtomSpec::new(Element::C)))
            .collect();
        let mut bonds: Vec<_> = (0..3)
            .map(|i| {
                Bond::from_spec(
                    BondId::new(i),
                    BondSpec::new(
                        AtomId::new(i),
                        AtomId::new(i + 1),
                        if i == 1 {
                            BondOrder::Double
                        } else {
                            BondOrder::Single
                        },
                    ),
                )
            })
            .collect();
        bonds[1].set_stereo_atoms(Some([AtomId::new(0), AtomId::new(3)]));
        bonds[1].set_stereo(BondStereo::Cis).unwrap();
        TopologyBlock::try_from_parts(atoms, bonds, vec![], vec![]).unwrap()
    }
    fn uint(value: &PropertyValue) -> Result<u32, TopologyEditError> {
        value
            .as_uint()
            .map_err(|_| TopologyEditError::SourceGroupParentCycle {
                group: SubstanceGroupId::new(999),
            })
    }
    #[test]
    fn source_batch_remove_bonds_empty_or_absent_mask_reads_no_state() {
        for pending in [None, Some(&[][..]), Some(&[false, false, false][..])] {
            let mut topology = topology();
            let before = topology.clone();
            let mut bookmarks = BTreeMap::from([(1, vec![])]);
            let before_marks = bookmarks.clone();
            batch_remove_bonds_source::<TopologyEditError>(
                &mut topology,
                pending,
                Some(&mut bookmarks),
                &mut |_| panic!("no source property read"),
            )
            .unwrap();
            assert_eq!(topology, before);
            assert_eq!(bookmarks, before_marks);
        }
    }
    #[test]
    fn source_batch_remove_bonds_descending_deletes_then_renumbers_aliases_and_groups() {
        let mut topology = topology();
        topology.stereo_groups = vec![
            StereoGroup::new(StereoGroupKind::And, vec![], vec![BondId::new(1)])
                .expect("valid distinct stereo members"),
        ];
        topology.substance_groups = vec![
            SubstanceGroup::new(SubstanceGroupId::new(0), SubstanceGroupKind::Data)
                .with_bonds(vec![BondId::new(1)]),
        ];
        let pending = [true, false, true];
        let mut bookmarks = BTreeMap::from([
            (1, vec![]),
            (2, vec![BondId::new(1)]),
            (3, vec![BondId::new(2)]),
        ]);
        batch_remove_bonds_source(
            &mut topology,
            Some(&pending),
            Some(&mut bookmarks),
            &mut uint,
        )
        .unwrap();
        assert_eq!(pending, [true, false, true]);
        assert_eq!(topology.bonds.len(), 1);
        assert_eq!(topology.bonds[0].id(), BondId::new(0));
        assert_eq!(topology.bonds[0].stereo(), BondStereo::None);
        assert_eq!(topology.bonds[0].stereo_atoms(), None);
        assert_eq!(topology.substance_groups[0].bonds(), [BondId::new(0)]);
        assert_eq!(topology.stereo_groups[0].bonds(), [BondId::new(0)]);
        assert_eq!(
            bookmarks,
            BTreeMap::from([(1, vec![]), (2, vec![BondId::new(0)])])
        );
        topology.validate().unwrap();
    }
    #[test]
    fn source_batch_remove_bonds_later_group_error_retains_prior_high_deletion_and_unrenumbered_ids()
     {
        let mut topology = topology();
        let mut group = SubstanceGroup::new(SubstanceGroupId::new(0), SubstanceGroupKind::Data)
            .with_bonds(vec![BondId::new(0)]);
        group.set_prop("PARENT", PropertyValue::UInt(10)).unwrap();
        group.set_prop("index", PropertyValue::Bool(false)).unwrap();
        topology.substance_groups.push(group);
        let pending = [true, false, true];
        let mut bookmarks = BTreeMap::from([
            (1, vec![BondId::new(0)]),
            (2, vec![BondId::new(2)]),
            (3, vec![BondId::new(1)]),
        ]);
        let error = batch_remove_bonds_source(
            &mut topology,
            Some(&pending),
            Some(&mut bookmarks),
            &mut uint,
        )
        .unwrap_err();
        assert!(matches!(
            error,
            TopologyEditError::SourceGroupParentCycle { .. }
        ));
        assert_eq!(topology.bonds.len(), 2); // high index2 deleted before low index0 failure
        assert_eq!(topology.bonds[0].id(), BondId::new(0));
        assert_eq!(topology.bonds[1].id(), BondId::new(1)); // final renumber not reached
        assert_eq!(topology.bonds[1].stereo(), BondStereo::None);
        assert_eq!(bookmarks, BTreeMap::from([(3, vec![BondId::new(1)])]));
        assert_eq!(topology.substance_groups[0].bonds(), [BondId::new(0)]);
    }
}

fn clear_source_atom_bookmarks(
    bookmarks: &mut std::collections::BTreeMap<i32, Vec<AtomId>>,
    atom: AtomId,
) {
    // RDKit❗✔️: if (std::find(atoms.begin(), atoms.end(), atom) != atoms.end()) {
    // RDKit❗✔️:   clearAtomBookmark(tmpI->first, atom);
    // RDKit❗✔️: }
    bookmarks.retain(|_, refs| {
        if let Some(position) = refs.iter().position(|id| *id == atom) {
            refs.remove(position);
            !refs.is_empty()
        } else {
            true
        }
    });
}

#[allow(dead_code)]
fn batch_remove_atoms_source<E: From<TopologyEditError>>(
    topology: &mut TopologyBlock,
    pending: Option<&[bool]>,
    mut bookmarks: Option<&mut std::collections::BTreeMap<i32, Vec<AtomId>>>,
    coordinates: &mut crate::CoordinateBlock,
    uint_reader: &mut (impl FnMut(&crate::PropertyValue) -> Result<u32, E> + ?Sized),
    text_reader: &mut (impl FnMut(&crate::PropertyValue) -> Result<crate::PropertyText, E> + ?Sized),
) -> Result<(), E> {
    // RDKit❗❌: void RWMol::batchRemoveAtoms() {
    // RDKit❗❌:   if (!dp_delAtoms || dp_delAtoms->none()) {
    // RDKit❗❌:     return;
    // RDKit❗❌:   }
    // RDKit❗❌:   std::vector<Atom *> oldIndices(getNumAtoms());
    // RDKit❗❌:   for (auto *atom : atoms()) {
    // RDKit❗❌:     oldIndices[atom->getIdx()] = atom;
    // RDKit❗❌:   }
    // RDKit❗❌:
    // RDKit❗❌:   auto &delAtoms = *dp_delAtoms;
    // RDKit❗❌:   for (unsigned int i = rdcast<unsigned int>(delAtoms.size()); i > 0; --i) {
    // RDKit❗❌:     if (!delAtoms[i - 1]) {
    // RDKit❗❌:       continue;
    // RDKit❗❌:     }
    // RDKit❗❌:     unsigned int idx = i - 1;
    // RDKit❗❌:     Atom *atom = getAtomWithIdx(idx);
    // RDKit❗❌:     if (!atom) {
    // RDKit❗❌:       continue;
    // RDKit❗❌:     }
    // RDKit❗❌:
    // RDKit❗❌:     // remove any bookmarks which point to this atom:
    // RDKit❗❌:     ATOM_BOOKMARK_MAP *marks = getAtomBookmarks();
    // RDKit❗❌:     auto markI = marks->begin();
    // RDKit❗❌:     while (markI != marks->end()) {
    // RDKit❗❌:       const ATOM_PTR_LIST &atoms = markI->second;
    // RDKit❗❌:       // we need to copy the iterator then increment it, because the
    // RDKit❗❌:       // deletion we're going to do in clearAtomBookmark will invalidate
    // RDKit❗❌:       // it.
    // RDKit❗❌:       auto tmpI = markI;
    // RDKit❗❌:       ++markI;
    // RDKit❗❌:       if (std::find(atoms.begin(), atoms.end(), atom) != atoms.end()) {
    // RDKit❗❌:         clearAtomBookmark(tmpI->first, atom);
    // RDKit❗❌:       }
    // RDKit❗❌:     }
    // RDKit❗❌:
    // RDKit❗❌:     // now deal with bonds:
    // RDKit❗❌:     //   their end indices may need to be decremented and their
    // RDKit❗❌:     //   indices will need to be handled and if they have an
    // RDKit❗❌:     //   ENDPTS prop that includes idx, it will need updating.
    // RDKit❗❌:     EDGE_ITER beg, end;
    // RDKit❗❌:     boost::tie(beg, end) = getEdges();
    // RDKit❗❌:     std::string sprop;
    // RDKit❗❌:     while (beg != end) {
    // RDKit❗❌:       Bond *bond = d_graph[*beg++];
    // RDKit❗❌:       if (bond->getPropIfPresent(RDKit::common_properties::_MolFileBondEndPts,
    // RDKit❗❌:                                  sprop)) {
    // RDKit❗❌:         // This would ideally use ParseV3000Array but I'm buggered if I can
    // RDKit❗❌:         // get the linker to find it.
    // RDKit❗❌:         //      std::vector<unsigned int> oats =
    // RDKit❗❌:         //          RDKit::SGroupParsing::ParseV3000Array<unsigned
    // RDKit❗❌:         //          int>(sprop);
    // RDKit❗❌:         if ('(' == sprop.front() && ')' == sprop.back()) {
    // RDKit❗❌:           sprop = sprop.substr(1, sprop.length() - 2);
    // RDKit❗❌:
    // RDKit❗❌:           // This is doing what ParseV3000Array would do.
    // RDKit❗❌:           boost::char_separator<char> sep(" ");
    // RDKit❗❌:           boost::tokenizer<boost::char_separator<char>> tokens(sprop, sep);
    // RDKit❗❌:           unsigned int num_ats = std::stod(*tokens.begin());
    // RDKit❗❌:           std::vector<unsigned int> oats;
    // RDKit❗❌:           auto beg = tokens.begin();
    // RDKit❗❌:           ++beg;
    // RDKit❗❌:           std::transform(beg, tokens.end(), std::back_inserter(oats),
    // RDKit❗❌:                          [](const std::string &a) { return std::stod(a); });
    // RDKit❗❌:
    // RDKit❗❌:           auto idx_pos = std::find(oats.begin(), oats.end(), idx + 1);
    // RDKit❗❌:           if (idx_pos != oats.end()) {
    // RDKit❗❌:             oats.erase(idx_pos);
    // RDKit❗❌:             --num_ats;
    // RDKit❗❌:           }
    // RDKit❗❌:           if (!num_ats) {
    // RDKit❗❌:             bond->clearProp(RDKit::common_properties::_MolFileBondEndPts);
    // RDKit❗❌:             bond->clearProp(common_properties::_MolFileBondAttach);
    // RDKit❗❌:           } else {
    // RDKit❗❌:             sprop = "(" + std::to_string(num_ats) + " ";
    // RDKit❗❌:             for (auto &i : oats) {
    // RDKit❗❌:               if (i > idx + 1) {
    // RDKit❗❌:                 --i;
    // RDKit❗❌:               }
    // RDKit❗❌:               sprop += std::to_string(i) + " ";
    // RDKit❗❌:             }
    // RDKit❗❌:             sprop[sprop.length() - 1] = ')';
    // RDKit❗❌:             bond->setProp(RDKit::common_properties::_MolFileBondEndPts, sprop);
    // RDKit❗❌:           }
    // RDKit❗❌:         }
    // RDKit❗❌:       }
    // RDKit❗❌:     }
    // RDKit❗❌:
    // RDKit❗❌:     removeSubstanceGroupsReferencingAtom(*this, idx);
    // RDKit❗❌:
    // RDKit❗❌:     // Remove this atom from any stereo group
    // RDKit❗❌:     removeAtomFromGroups(atom, d_stereo_groups);
    // RDKit❗❌:     atom->setOwningMol(nullptr);
    // RDKit❗❌:
    // RDKit❗❌:     // remove all connections to the atom:
    // RDKit❗❌:     MolGraph::vertex_descriptor vd = boost::vertex(idx, d_graph);
    // RDKit❗❌:     boost::clear_vertex(vd, d_graph);
    // RDKit❗❌:     // finally remove the vertex itself
    // RDKit❗❌:     boost::remove_vertex(vd, d_graph);
    // RDKit❗❌:     delete atom;
    // RDKit❗❌:     oldIndices[idx] = nullptr;
    // RDKit❗❌:   }
    // RDKit❗❌:
    // RDKit❗❌:   // reassign atom indices
    // RDKit❗❌:   for (unsigned int i = 0; i < getNumAtoms(); i++) {
    // RDKit❗❌:     Atom *atm = getAtomWithIdx(i);
    // RDKit❗❌:     atm->setIdx(i);
    // RDKit❗❌:   }
    // RDKit❗❌:
    // RDKit❗❌:   // reassign atom indices in bonds
    // RDKit❗❌:   for (auto bond : bonds()) {
    // RDKit❗❌:     auto bgnidx = bond->getBeginAtomIdx();
    // RDKit❗❌:     auto endidx = bond->getEndAtomIdx();
    // RDKit❗❌:     Atom *bgn = oldIndices[bgnidx];
    // RDKit❗❌:     Atom *end = oldIndices[endidx];
    // RDKit❗❌:     CHECK_INVARIANT(bgn, "Atom mapping failed");
    // RDKit❗❌:     CHECK_INVARIANT(end, "Atom mapping failed");
    // RDKit❗❌:     bond->setBeginAtomIdx(bgn->getIdx());
    // RDKit❗❌:     bond->setEndAtomIdx(end->getIdx());
    // RDKit❗❌:     INT_VECT stereoAtoms;
    // RDKit❗❌:     INT_VECT &oldStereoAtoms = bond->getStereoAtoms();
    // RDKit❗❌:     if (oldStereoAtoms.size()) {
    // RDKit❗❌:       for (auto &idx : oldStereoAtoms) {
    // RDKit❗❌:         if (oldIndices[idx]) {
    // RDKit❗❌:           stereoAtoms.push_back(oldIndices[idx]->getIdx());
    // RDKit❗❌:         }
    // RDKit❗❌:       }
    // RDKit❗❌:       bond->getStereoAtoms().swap(stereoAtoms);
    // RDKit❗❌:     }
    // RDKit❗❌:   }
    // RDKit❗❌:
    // RDKit❗❌:   // do the same with the coordinates in the conformations
    // RDKit❗❌:   for (auto conf : d_confs) {
    // RDKit❗❌:     RDGeom::POINT3D_VECT &positions = conf->getPositions();
    // RDKit❗❌:     RDGeom::POINT3D_VECT newPositions;
    // RDKit❗❌:     newPositions.reserve(getNumAtoms());
    // RDKit❗❌:
    // RDKit❗❌:     for (RDGeom::POINT3D_VECT::size_type i = 0; i < positions.size(); ++i) {
    // RDKit❗❌:       if (oldIndices[i] != nullptr) {
    // RDKit❗❌:         newPositions.push_back(positions[i]);
    // RDKit❗❌:       }
    // RDKit❗❌:     }
    // RDKit❗❌:     CHECK_INVARIANT(newPositions.size() == getNumAtoms(), "Lost coordinates!");
    // RDKit❗❌:     positions.swap(newPositions);
    // RDKit❗❌:   }
    // RDKit❗❌: }
    // Native pointer/null/signed-width/allocator/dangling states are unmodeled.
    // This source-state kernel does not claim the structural mapping finish as
    // native batch commit. No ring/property resets occur in this source method.
    // Cost: Vec vertex/edge erasure, extra stable-ID aliases and CSR rebuilding
    // add source-absent work. ENDPTS/SG effects precede each physical deletion;
    // endpoint/reference and conformer remapping happen only after all deletions.
    let Some(pending) = pending else {
        return Ok(());
    };
    if !pending.iter().any(|removed| *removed) {
        return Ok(());
    }
    let original_count = topology.atoms.len();
    let mut aliases = vec![None; original_count];
    for index in (0..pending.len()).rev() {
        if !pending[index] {
            continue;
        }
        let atom = topology
            .atoms
            .get(index)
            .ok_or_else(|| {
                E::from(TopologyEditError::AtomOutOfRange {
                    atom: AtomId::new(index),
                    atom_count: topology.atoms.len(),
                })
            })?
            .id();
        if let Some(bookmarks) = bookmarks.as_deref_mut() {
            clear_source_atom_bookmarks(bookmarks, atom);
        }
        for bond in &mut topology.bonds {
            if let Some(value) = bond.prop("_MolFileBondEndPts") {
                let text = text_reader(value)?;
                update_source_bond_endpoints_property(bond, index, &text).map_err(E::from)?;
            }
        }
        crate::sgroup::remove_source_groups_referencing_atom(
            &mut topology.substance_groups,
            AtomId::new(index),
            uint_reader,
        )?;
        for group in &mut topology.stereo_groups {
            group.remove_atom(atom);
        }
        topology.stereo_groups.retain(|group| !group.is_empty());
        // clear_vertex removes all actual remaining graph connections. The
        // native standalone path can leave its separate bond counter/IDs stale;
        // derived model counts have no such independent counter state.
        topology
            .bonds
            .retain(|bond| bond.begin() != atom && bond.end() != atom);
        topology.atoms.remove(index);
    }
    for (index, atom) in topology.atoms.iter_mut().enumerate() {
        let old = atom.id();
        let target = aliases.get_mut(old.index()).ok_or_else(|| {
            E::from(TopologyEditError::AtomOutOfRange {
                atom: old,
                atom_count: original_count,
            })
        })?;
        atom.set_source_index(AtomId::new(index));
        *target = Some(atom.id());
    }
    // Retained native pointers see new atom IDs at setIdx, including group and
    // bookmark aliases. First-only duplicate deleted aliases remain dangling
    // native state; they are not silently filtered or rebound here.
    if let Some(bookmarks) = bookmarks {
        for values in bookmarks.values_mut() {
            for id in values {
                if let Some(Some(mapped)) = aliases.get(id.index()) {
                    *id = *mapped;
                }
            }
        }
    }
    for group in &mut topology.stereo_groups {
        group.remap_source_retained_atom_ids(&aliases);
    }
    for bond in &mut topology.bonds {
        let resolve = |endpoint: &'static str, atom: AtomId| -> Result<AtomId, E> {
            aliases.get(atom.index()).copied().flatten().ok_or_else(|| {
                E::from(TopologyEditError::InvalidSource(
                    TopologyValidationError::BondEndpointOutOfRange {
                        bond: bond.id(),
                        endpoint,
                        atom,
                        atom_count: original_count,
                    },
                ))
            })
        };
        let begin = resolve("begin", bond.begin())?;
        let end = resolve("end", bond.end())?;
        bond.set_endpoints(begin, end);
        if !bond.stereo_atom_references().is_empty() {
            let mut references = Vec::new();
            for atom in bond.stereo_atom_references() {
                let mapped = aliases.get(atom.index()).ok_or_else(|| {
                    E::from(TopologyEditError::InvalidSource(
                        TopologyValidationError::StereoReferenceOutOfRange {
                            bond: bond.id(),
                            atom: *atom,
                            atom_count: original_count,
                        },
                    ))
                })?;
                if let Some(mapped) = mapped {
                    references.push(*mapped);
                }
            }
            // Native filters each entry independently: one surviving reference
            // stays a one-element vector; stereo code itself is not cleared.
            bond.set_source_stereo_atom_references(references);
        }
    }
    topology.adjacency = AdjacencyList::from_topology(topology.atoms.len(), &topology.bonds);
    coordinates
        .source_batch_remove_atom_positions(&aliases, topology.atoms.len())
        .map_err(|error| E::from(TopologyEditError::Coordinate(error)))?;
    Ok(())
}

#[cfg(test)]
mod source_batch_remove_atoms_complete_tests {
    use super::*;
    use crate::{
        BondOrder, Conformer2D, Conformer3D, CoordinateBlock, CoordinateDimension, Element,
        PropertyValue, StereoGroupKind, SubstanceGroupId, SubstanceGroupKind,
    };
    use std::collections::BTreeMap;
    fn topology() -> TopologyBlock {
        let atoms = (0..4)
            .map(|i| Atom::from_spec(AtomId::new(i), AtomSpec::new(Element::C)))
            .collect();
        let mut bond = Bond::from_spec(
            BondId::new(0),
            BondSpec::new(AtomId::new(0), AtomId::new(2), BondOrder::Double)
                .with_stereo_atoms(AtomId::new(1), AtomId::new(3)),
        );
        bond.set_stereo(BondStereo::Cis).unwrap();
        TopologyBlock::try_from_parts(atoms, vec![bond], vec![], vec![]).unwrap()
    }
    fn uint(value: &PropertyValue) -> Result<u32, TopologyEditError> {
        value
            .as_uint()
            .map_err(|_| TopologyEditError::SourceGroupParentCycle {
                group: SubstanceGroupId::new(999),
            })
    }
    fn text(value: &PropertyValue) -> Result<crate::PropertyText, TopologyEditError> {
        value
            .as_string()
            .cloned()
            .map_err(|_| TopologyEditError::SourceGroupParentCycle {
                group: SubstanceGroupId::new(998),
            })
    }
    #[test]
    fn source_batch_remove_atoms_absent_or_empty_mask_reads_no_state() {
        for mask in [None, Some(&[][..]), Some(&[false; 4][..])] {
            let mut graph = topology();
            let before = graph.clone();
            let mut coords = CoordinateBlock::default();
            let mut marks = BTreeMap::from([(1, vec![])]);
            batch_remove_atoms_source::<TopologyEditError>(
                &mut graph,
                mask,
                Some(&mut marks),
                &mut coords,
                &mut |_| panic!("no property read"),
                &mut |_| panic!("no property read"),
            )
            .unwrap();
            assert_eq!(graph, before);
            assert_eq!(marks, BTreeMap::from([(1, vec![])]));
        }
    }
    #[test]
    fn source_batch_remove_atoms_preserves_single_reference_and_stereo_tag_with_final_coordinates()
    {
        let mut graph = topology();
        graph.bonds[0]
            .set_prop(
                "_MolFileBondEndPts",
                PropertyValue::String("(3 2 4 4)".into()),
            )
            .unwrap();
        graph.stereo_groups = vec![
            StereoGroup::new(StereoGroupKind::And, vec![AtomId::new(3)], vec![])
                .expect("valid distinct stereo members"),
        ];
        graph.substance_groups = vec![
            SubstanceGroup::new(SubstanceGroupId::new(0), SubstanceGroupKind::Data)
                .with_atoms(vec![AtomId::new(3)]),
        ];
        let negative_zero = -0.0f64;
        let mut coords = CoordinateBlock {
            conformers_2d: vec![
                Conformer2D::new(
                    4,
                    vec![[negative_zero, 1.0], [2.0, 3.0], [4.0, 5.0], [6.0, 7.0]],
                )
                .with_prop("note", "keep"),
            ],
            conformers_3d: vec![Conformer3D::new(
                8,
                vec![[0.0; 3], [1.0; 3], [2.0; 3], [3.0; 3]],
                false,
            )],
            source_conformer_order: Some(vec![
                CoordinateDimension::ThreeD,
                CoordinateDimension::TwoD,
            ]),
            ..Default::default()
        };
        let before_props = coords.conformers_2d[0].props().clone();
        let mut marks = BTreeMap::from([
            (1, vec![]),
            (2, vec![AtomId::new(1)]),
            (3, vec![AtomId::new(3)]),
        ]);
        batch_remove_atoms_source(
            &mut graph,
            Some(&[false, true, false, false]),
            Some(&mut marks),
            &mut coords,
            &mut uint,
            &mut text,
        )
        .unwrap();
        assert_eq!(
            graph.atoms.iter().map(Atom::id).collect::<Vec<_>>(),
            vec![AtomId::new(0), AtomId::new(1), AtomId::new(2)]
        );
        assert_eq!(
            (graph.bonds[0].begin(), graph.bonds[0].end()),
            (AtomId::new(0), AtomId::new(1))
        );
        assert_eq!(graph.bonds[0].stereo_atom_references(), [AtomId::new(2)]);
        assert_eq!(graph.bonds[0].stereo(), BondStereo::Cis); // native does not clear this code
        assert_eq!(
            graph.bonds[0].prop("_MolFileBondEndPts"),
            Some(&PropertyValue::String("(2 3 3)".into()))
        );
        assert_eq!(graph.substance_groups[0].atoms(), [AtomId::new(2)]);
        assert_eq!(graph.stereo_groups[0].atoms(), [AtomId::new(2)]);
        assert_eq!(
            marks,
            BTreeMap::from([(1, vec![]), (3, vec![AtomId::new(2)])])
        );
        assert_eq!(
            coords.conformers_2d[0].coordinates(),
            [[negative_zero, 1.0], [4.0, 5.0], [6.0, 7.0]]
        );
        assert_eq!(
            coords.conformers_2d[0].coordinates()[0][0].to_bits(),
            negative_zero.to_bits()
        );
        assert_eq!(coords.conformers_2d[0].props(), &before_props);
        assert_eq!(
            coords.conformers_3d[0].coordinates(),
            [[0.0; 3], [2.0; 3], [3.0; 3]]
        );
        assert!(!coords.conformers_3d[0].is_3d());
        assert!(matches!(
            graph.validate(),
            Err(TopologyValidationError::StereoAtomsRequired { .. })
        )); // source-produced raw state is retained; local validation still rejects malformed CIS pair
    }
    #[test]
    fn source_batch_remove_atoms_filters_all_vector_cardinalities_in_order() {
        let mut graph = topology();
        graph.bonds[0].set_source_stereo_atom_references(vec![
            AtomId::new(3),
            AtomId::new(1),
            AtomId::new(3),
            AtomId::new(0),
        ]);
        let mut marks = BTreeMap::new();
        batch_remove_atoms_source(
            &mut graph,
            Some(&[false, true, false, false]),
            Some(&mut marks),
            &mut CoordinateBlock::default(),
            &mut uint,
            &mut text,
        )
        .unwrap();
        assert_eq!(
            graph.bonds[0].stereo_atom_references(),
            [AtomId::new(2), AtomId::new(2), AtomId::new(0)]
        );
        let mut clone = graph.bonds[0].clone();
        clone.set_source_stereo_atom_references(vec![AtomId::new(9)]);
        graph.bonds[0] = clone;
        assert!(
            matches!(graph.validate(), Err(TopologyValidationError::StereoReferenceOutOfRange { atom, .. }) if atom == AtomId::new(9))
        );
    }
    #[test]
    fn source_batch_remove_atoms_late_coordinate_error_retains_prior_frame_swap_and_graph_deletion()
    {
        let mut graph = topology();
        let mut marks = BTreeMap::new();
        let mut coords = CoordinateBlock {
            conformers_2d: vec![Conformer2D::new(
                1,
                vec![[0.0; 2], [1.0; 2], [2.0; 2], [3.0; 2]],
            )],
            conformers_3d: vec![Conformer3D::new(2, vec![[0.0; 3]], true)],
            source_conformer_order: Some(vec![
                CoordinateDimension::TwoD,
                CoordinateDimension::ThreeD,
            ]),
            ..Default::default()
        };
        let error = batch_remove_atoms_source(
            &mut graph,
            Some(&[false, true, false, false]),
            Some(&mut marks),
            &mut coords,
            &mut uint,
            &mut text,
        )
        .unwrap_err();
        assert!(matches!(
            error,
            TopologyEditError::Coordinate(crate::CoordinateValidationError::RowCount {
                dimension: "3D",
                rows: 1,
                atom_count: 3,
                ..
            })
        ));
        assert_eq!(graph.atoms.len(), 3);
        assert_eq!(graph.bonds[0].stereo_atom_references(), [AtomId::new(2)]);
        assert_eq!(
            coords.conformers_2d[0].coordinates(),
            [[0.0; 2], [2.0; 2], [3.0; 2]]
        );
        assert_eq!(coords.conformers_3d[0].coordinates(), [[0.0; 3]]);
    }
    #[test]
    fn source_batch_remove_atoms_accepts_short_coordinates_when_only_deleted_tail_is_missing() {
        let mut graph = topology();
        let mut marks = BTreeMap::new();
        let mut coords = CoordinateBlock {
            conformers_2d: vec![Conformer2D::new(3, vec![[0.0; 2], [1.0; 2], [2.0; 2]])],
            ..Default::default()
        };
        batch_remove_atoms_source(
            &mut graph,
            Some(&[false, false, false, true]),
            Some(&mut marks),
            &mut coords,
            &mut uint,
            &mut text,
        )
        .unwrap();
        assert_eq!(coords.conformers_2d[0].coordinates().len(), 3);
        assert_eq!(graph.bonds[0].stereo_atom_references(), [AtomId::new(1)]);
    }
    #[test]
    fn source_batch_remove_atoms_later_property_error_retains_descending_deletion_before_final_remapping()
     {
        let mut graph = topology();
        graph.bonds[0].set_endpoints(AtomId::new(0), AtomId::new(3));
        graph.bonds[0]
            .set_prop("_MolFileBondEndPts", PropertyValue::String("(1 1)".into()))
            .unwrap();
        let mut marks = BTreeMap::from([
            (1, vec![AtomId::new(1)]),
            (2, vec![AtomId::new(2)]),
            (3, vec![AtomId::new(3)]),
        ]);
        let mut coords = CoordinateBlock {
            conformers_2d: vec![Conformer2D::new(
                1,
                vec![[0.0; 2], [1.0; 2], [2.0; 2], [3.0; 2]],
            )],
            ..Default::default()
        };
        let before_coords = coords.clone();
        let mut calls = 0;
        let error = batch_remove_atoms_source(
            &mut graph,
            Some(&[false, true, true, false]),
            Some(&mut marks),
            &mut coords,
            &mut uint,
            &mut |value| {
                calls += 1;
                if calls == 2 {
                    return Err(TopologyEditError::SourceGroupParentCycle {
                        group: SubstanceGroupId::new(998),
                    });
                }
                text(value)
            },
        )
        .unwrap_err();
        assert_eq!(calls, 2);
        assert!(
            matches!(error, TopologyEditError::SourceGroupParentCycle { group } if group == SubstanceGroupId::new(998))
        );
        assert_eq!(
            graph.atoms.iter().map(Atom::id).collect::<Vec<_>>(),
            vec![AtomId::new(0), AtomId::new(1), AtomId::new(3)]
        );
        assert_eq!(graph.bonds[0].end(), AtomId::new(3));
        assert_eq!(
            graph.bonds[0].stereo_atom_references(),
            [AtomId::new(1), AtomId::new(3)]
        );
        assert_eq!(marks, BTreeMap::from([(3, vec![AtomId::new(3)])]));
        assert_eq!(coords, before_coords); // no coordinates changed until all deletes complete
    }
}

fn clear_source_computed_properties<E: From<TopologyEditError>>(
    topology: &mut TopologyBlock,
    properties: &mut crate::MoleculeProperties,
    include_rings: bool,
    reset_ring: Option<&mut dyn FnMut()>,
) -> Result<(), E> {
    // RDKit❗❌: void ROMol::clearComputedProps(bool includeRings) const {
    // RDKit❗❌:   // the SSSR information:
    // RDKit❗❌:   if (includeRings) {
    // RDKit❗❌:     this->dp_ringInfo->reset();
    // RDKit❗❌:   }
    // RDKit❗❌:
    // RDKit❗❌:   RDProps::clearComputedProps();
    // RDKit❗❌:
    // RDKit❗❌:   for (auto atom : atoms()) {
    // RDKit❗❌:     atom->clearComputedProps();
    // RDKit❗❌:   }
    // RDKit❗❌:
    // RDKit❗❌:   for (auto bond : bonds()) {
    // RDKit❗❌:     bond->clearComputedProps();
    // RDKit❗❌:   }
    // RDKit❗❌: }
    // One canonical molecule/atom/bond computed-property operation; preserve
    // source order and every earlier mutation on a later property-kind error.
    // Native RingInfo object identity/ownership is supplied through this actual
    // state callback, not an invented dummy ring value. Typed dictionary costs
    // and source computed-list conversion gaps remain separately recorded.
    if include_rings && let Some(reset_ring) = reset_ring {
        reset_ring();
    }
    properties
        .clear_computed_props()
        .map_err(|error| E::from(TopologyEditError::MoleculeProperty(error)))?;
    for atom in &mut topology.atoms {
        atom.clear_computed_props()
            .map_err(|error| E::from(TopologyEditError::AtomProperty(error)))?;
    }
    for bond in &mut topology.bonds {
        bond.clear_computed_props()
            .map_err(|error| E::from(TopologyEditError::InvalidBond(error)))?;
    }
    Ok(())
}

// Explicit borrowed detached source state: masks/bookmarks/coordinates are the
// caller's actual values; this view stores no second graph or live authority.
/// Borrowed actual detached source state for native batch deletion.
/// `None` metadata is independently unmodeled, not a source empty/default value;
/// full native bookmark/cache capability requires actual `Some` values.
#[doc(hidden)]
pub struct SourceBatchCommitState<'a, E> {
    pub atoms: &'a mut Option<Vec<bool>>,
    pub bonds: &'a mut Option<Vec<bool>>,
    pub atom_bookmarks: Option<&'a mut std::collections::BTreeMap<i32, Vec<AtomId>>>,
    pub bond_bookmarks: Option<&'a mut std::collections::BTreeMap<i32, Vec<BondId>>>,
    pub coordinates: &'a mut crate::CoordinateBlock,
    pub properties: &'a mut crate::MoleculeProperties,
    pub reset_ring: Option<&'a mut dyn FnMut()>,
    pub uint_reader: &'a mut dyn FnMut(&crate::PropertyValue) -> Result<u32, E>,
    pub text_reader: &'a mut dyn FnMut(&crate::PropertyValue) -> Result<crate::PropertyText, E>,
}

#[allow(dead_code)]
pub fn commit_batch_edit_source<E: From<TopologyEditError>>(
    topology: &mut TopologyBlock,
    state: SourceBatchCommitState<'_, E>,
) -> Result<(), E> {
    // RDKit❗❌: void RWMol::commitBatchEdit() {
    // RDKit❗❌:   if (!(dp_delBonds || dp_delAtoms)) {
    // RDKit❗❌:     return;
    // RDKit❗❌:   } else if (dp_delBonds->none() && dp_delAtoms->none()) {
    // RDKit❗❌:     // no need to reset ring info & calculated properties,
    // RDKit❗❌:     // since nothing gets removed
    // RDKit❗❌:     dp_delBonds.reset();
    // RDKit❗❌:     dp_delAtoms.reset();
    // RDKit❗❌:     return;
    // RDKit❗❌:   }
    // RDKit❗❌:
    // RDKit❗❌:   batchRemoveBonds();
    // RDKit❗❌:   batchRemoveAtoms();
    // RDKit❗❌:
    // RDKit❗❌:   // remove ring info
    // RDKit❗❌:   dp_ringInfo->reset();
    // RDKit❗❌:
    // RDKit❗❌:   // fix properties
    // RDKit❗❌:   clearComputedProps(true);
    // RDKit❗❌:   dp_delBonds.reset();
    // RDKit❗❌:   dp_delAtoms.reset();
    // RDKit❗❌: }
    // Source commit is separate from local structural mapping finish. Reuse
    // the complete canonical batch methods rather than bulk-filtering edges,
    // atoms and all annotations simultaneously. Native nullable partial mask
    // state would dereference null; translate it to a structural state error.
    // This private source-state entry is integrated/reviewed at the whole gate.
    let SourceBatchCommitState {
        atoms,
        bonds,
        atom_bookmarks,
        bond_bookmarks,
        coordinates,
        properties,
        mut reset_ring,
        uint_reader,
        text_reader,
    } = state;
    match (bonds.as_deref(), atoms.as_deref()) {
        (None, None) => return Ok(()),
        (Some(bonds_mask), Some(atoms_mask)) => {
            if !bonds_mask.iter().any(|flag| *flag) && !atoms_mask.iter().any(|flag| *flag) {
                *bonds = None;
                *atoms = None;
                return Ok(());
            }
        }
        _ => return Err(TopologyEditError::IncompleteSourceBatchMasks.into()),
    }
    batch_remove_bonds_source(topology, bonds.as_deref(), bond_bookmarks, uint_reader)?;
    batch_remove_atoms_source(
        topology,
        atoms.as_deref(),
        atom_bookmarks,
        coordinates,
        uint_reader,
        text_reader,
    )?;
    if let Some(reset_ring) = reset_ring.as_deref_mut() {
        reset_ring();
    }
    // Native explicitly resets once here, and clearComputedProps(true) resets
    // again. Do not deduplicate these observable calls or reset masks on error.
    clear_source_computed_properties(topology, properties, true, reset_ring)?;
    *bonds = None;
    *atoms = None;
    Ok(())
}

#[cfg(test)]
mod source_commit_batch_edit_complete_tests {
    use super::*;
    use crate::{
        BondOrder, Conformer2D, CoordinateBlock, Element, MoleculeProperties, PropertyValue,
    };
    use std::{cell::Cell, collections::BTreeMap};
    struct Fixture {
        graph: TopologyBlock,
        atoms: Option<Vec<bool>>,
        bonds: Option<Vec<bool>>,
        atom_marks: BTreeMap<i32, Vec<AtomId>>,
        bond_marks: BTreeMap<i32, Vec<BondId>>,
        coords: CoordinateBlock,
        props: MoleculeProperties,
        resets: Cell<usize>,
    }
    impl Fixture {
        fn new() -> Self {
            let atoms = (0..3)
                .map(|i| Atom::from_spec(AtomId::new(i), AtomSpec::new(Element::C)))
                .collect();
            let bonds = vec![
                Bond::from_spec(
                    BondId::new(0),
                    BondSpec::new(AtomId::new(0), AtomId::new(1), BondOrder::Single),
                ),
                Bond::from_spec(
                    BondId::new(1),
                    BondSpec::new(AtomId::new(0), AtomId::new(2), BondOrder::Single),
                ),
            ];
            let mut graph = TopologyBlock::try_from_parts(atoms, bonds, vec![], vec![]).unwrap();
            graph.atoms[0]
                .set_computed_prop("atom", PropertyValue::UInt(1))
                .unwrap();
            graph.bonds[1]
                .set_computed_prop("bond", PropertyValue::UInt(1))
                .unwrap();
            let mut props = MoleculeProperties::default();
            props
                .set_computed_prop("mol", PropertyValue::UInt(1))
                .unwrap();
            Self {
                graph,
                atoms: Some(vec![false, true, false]),
                bonds: Some(vec![true, false]),
                atom_marks: BTreeMap::from([(1, vec![AtomId::new(2)])]),
                bond_marks: BTreeMap::from([(1, vec![BondId::new(1)])]),
                coords: CoordinateBlock {
                    conformers_2d: vec![Conformer2D::new(3, vec![[0.0; 2], [1.0; 2], [2.0; 2]])],
                    ..Default::default()
                },
                props,
                resets: Cell::new(0),
            }
        }
        fn commit(&mut self) -> Result<(), TopologyEditError> {
            commit_batch_edit_source(
                &mut self.graph,
                SourceBatchCommitState {
                    atoms: &mut self.atoms,
                    bonds: &mut self.bonds,
                    atom_bookmarks: Some(&mut self.atom_marks),
                    bond_bookmarks: Some(&mut self.bond_marks),
                    coordinates: &mut self.coords,
                    properties: &mut self.props,
                    reset_ring: Some(&mut || self.resets.set(self.resets.get() + 1)),
                    uint_reader: &mut |_| panic!("fixture contains no SG properties"),
                    text_reader: &mut |_| panic!("fixture contains no ENDPTS properties"),
                },
            )
        }
    }
    #[test]
    fn source_commit_batch_edit_inactive_and_empty_masks_do_not_touch_ring_or_properties() {
        for inactive in [true, false] {
            let mut f = Fixture::new();
            let before = f.graph.clone();
            let coords = f.coords.clone();
            let props = f.props.clone();
            if inactive {
                f.atoms = None;
                f.bonds = None;
            } else {
                f.atoms = Some(vec![false; 3]);
                f.bonds = Some(vec![false; 2]);
            }
            f.commit().unwrap();
            assert_eq!(f.resets.get(), 0);
            assert_eq!(f.atoms, None);
            assert_eq!(f.bonds, None);
            assert_eq!(f.graph, before);
            assert_eq!(f.coords, coords);
            assert_eq!(f.props, props);
        }
    }
    #[test]
    fn source_commit_batch_edit_deletes_bonds_before_atoms_resets_twice_and_clears_properties() {
        let mut f = Fixture::new();
        f.graph.atoms[1]
            .set_prop("__computedProps", PropertyValue::Int(7))
            .unwrap(); // removed atom is never visited by final clear
        f.commit().unwrap();
        assert_eq!(f.resets.get(), 2);
        assert_eq!(f.atoms, None);
        assert_eq!(f.bonds, None);
        assert_eq!(f.graph.atoms.len(), 2);
        assert_eq!(f.graph.bonds.len(), 1);
        assert_eq!(
            (
                f.graph.bonds[0].id(),
                f.graph.bonds[0].begin(),
                f.graph.bonds[0].end()
            ),
            (BondId::new(0), AtomId::new(0), AtomId::new(1))
        );
        assert_eq!(f.atom_marks[&1], [AtomId::new(1)]);
        assert_eq!(f.bond_marks[&1], [BondId::new(0)]);
        assert_eq!(
            f.coords.conformers_2d[0].coordinates(),
            [[0.0; 2], [2.0; 2]]
        );
        assert!(f.props.prop("mol").is_none());
        assert!(f.graph.atoms[0].prop("atom").is_none());
        assert!(f.graph.bonds[0].prop("bond").is_none());
        f.graph.validate().unwrap();
        f.commit().unwrap();
        assert_eq!(f.resets.get(), 2); // successful masks are inactive
    }
    #[test]
    fn source_commit_batch_edit_coordinate_error_keeps_masks_and_skips_ring_and_property_clear() {
        let mut f = Fixture::new();
        f.coords.conformers_2d[0] = Conformer2D::new(3, vec![[0.0; 2]]);
        let error = f.commit().unwrap_err();
        assert!(matches!(
            error,
            TopologyEditError::Coordinate(crate::CoordinateValidationError::RowCount {
                rows: 1,
                atom_count: 2,
                ..
            })
        ));
        assert_eq!(f.graph.atoms.len(), 2);
        assert_eq!(f.graph.bonds.len(), 1);
        assert_eq!(f.resets.get(), 0);
        assert_eq!(f.atoms, Some(vec![false, true, false]));
        assert_eq!(f.bonds, Some(vec![true, false]));
        assert!(f.props.prop("mol").is_some());
        assert!(f.graph.atoms[0].prop("atom").is_some());
    }
    #[test]
    fn source_commit_batch_edit_property_error_keeps_masks_after_both_resets_and_prior_clears() {
        let mut f = Fixture::new();
        f.graph.bonds[1]
            .set_prop("__computedProps", PropertyValue::Int(7))
            .unwrap();
        let error = f.commit().unwrap_err();
        assert!(matches!(
            error,
            TopologyEditError::InvalidBond(BondValueError::ComputedListKind(_))
        ));
        assert_eq!(f.resets.get(), 2);
        assert_eq!(f.atoms, Some(vec![false, true, false]));
        assert_eq!(f.bonds, Some(vec![true, false]));
        assert!(f.props.prop("mol").is_none());
        assert!(f.graph.atoms[0].prop("atom").is_none());
        assert!(f.graph.bonds[0].prop("bond").is_some());
        assert_eq!(f.coords.conformers_2d[0].coordinates().len(), 2);
    }
    #[test]
    fn source_commit_batch_edit_partial_mask_state_reports_structural_error_before_mutation() {
        let mut f = Fixture::new();
        f.bonds = None;
        let before = f.graph.clone();
        assert_eq!(
            f.commit(),
            Err(TopologyEditError::IncompleteSourceBatchMasks)
        );
        assert_eq!(f.graph, before);
        assert_eq!(f.resets.get(), 0);
        assert_eq!(f.atoms, Some(vec![false, true, false]));
    }
}

/// Canonical source replacement of one detached bond, without a batch clone.
#[doc(hidden)]
pub fn replace_source_bond<E: From<TopologyEditError>>(
    topology: &mut TopologyBlock,
    bond: BondId,
    replacement: &Bond,
    preserve_props: bool,
    keep_sgroups: bool,
    mut order_reader: impl FnMut(crate::BondOrder) -> Result<f64, E>,
    mut uint_reader: impl FnMut(&crate::PropertyValue) -> Result<u32, E>,
) -> Result<(), E> {
    // RDKit❗❌: void RWMol::replaceBond(unsigned int idx, Bond *bond_pin, bool preserveProps,
    // RDKit❗❌:                         bool keepSGroups) {
    // RDKit❗❌:   PRECONDITION(bond_pin, "bad bond passed to replaceBond");
    // RDKit❗❌:   URANGE_CHECK(idx, getNumBonds());
    // RDKit❗❌:   auto bIter = getEdges();
    // RDKit❗❌:   for (unsigned int i = 0; i < idx; i++) {
    // RDKit❗❌:     ++bIter.first;
    // RDKit❗❌:   }
    // RDKit❗❌:   const auto *obond = d_graph[*(bIter.first)];
    // RDKit❗❌:   auto *bond_p = bond_pin->copy();
    // RDKit❗❌:   bond_p->setOwningMol(this);
    // RDKit❗❌:   bond_p->setIdx(idx);
    // RDKit❗❌:   bond_p->setBeginAtomIdx(obond->getBeginAtomIdx());
    // RDKit❗❌:   bond_p->setEndAtomIdx(obond->getEndAtomIdx());
    // RDKit❗❌:
    // RDKit❗❌:   // Update explicit Hs, if set, on both ends. This was github #7128
    // RDKit❗❌:   auto orderDifference =
    // RDKit❗❌:       bond_p->getBondTypeAsDouble() - obond->getBondTypeAsDouble();
    // RDKit❗❌:   if (orderDifference > 0) {
    // RDKit❗❌:     for (auto atom : {bond_p->getBeginAtom(), bond_p->getEndAtom()}) {
    // RDKit❗❌:       if (auto explicit_hs = atom->getNumExplicitHs(); explicit_hs > 0) {
    // RDKit❗❌:         auto new_hs = static_cast<int>(explicit_hs - orderDifference);
    // RDKit❗❌:         atom->setNumExplicitHs(std::max(new_hs, 0));
    // RDKit❗❌:       }
    // RDKit❗❌:     }
    // RDKit❗❌:   }
    // RDKit❗❌:
    // RDKit❗❌:   if (preserveProps) {
    // RDKit❗❌:     const bool replaceExistingData = false;
    // RDKit❗❌:     bond_p->updateProps(*d_graph[*(bIter.first)], replaceExistingData);
    // RDKit❗❌:   }
    // RDKit❗❌:
    // RDKit❗❌:   const auto orig_p = d_graph[*(bIter.first)];
    // RDKit❗❌:   delete orig_p;
    // RDKit❗❌:   d_graph[*(bIter.first)] = bond_p;
    // RDKit❗❌:
    // RDKit❗❌:   if (!keepSGroups) {
    // RDKit❗❌:     removeSubstanceGroupsReferencingBond(*this, idx);
    // RDKit❗❌:   }
    // RDKit❗❌:
    // RDKit❗❌:   // handle bookmarks
    // RDKit❗❌:   for (auto &ab : d_bondBookmarks) {
    // RDKit❗❌:     for (auto &elem : ab.second) {
    // RDKit❗❌:       if (elem == orig_p) {
    // RDKit❗❌:         elem = bond_p;
    // RDKit❗❌:       }
    // RDKit❗❌:     }
    // RDKit❗❌:   }
    // RDKit❗❌: }
    // Behavior gaps: explicit-H u8 versus native uint32; stable IDs versus
    // owning/bookmark pointers; native moved-from group state after errors.
    // Cost: indexed edge access is cheaper, but canonical property tree
    // cloning and detached SGroup identity remapping add work.
    let bond_count = topology.bonds.len();
    let old = topology
        .bonds
        .get(bond.index())
        .ok_or_else(|| E::from(TopologyEditError::BondOutOfRange { bond, bond_count }))?;
    let mut copied = replacement.clone();
    copied.set_id_for_construction(bond);
    copied.set_endpoints(old.begin(), old.end());
    // Preserve evaluation order: new type first, then original type.
    let order_difference = order_reader(copied.order())? - order_reader(old.order())?;
    if order_difference > 0.0 {
        for endpoint in [copied.begin(), copied.end()] {
            let atom = &mut topology.atoms[endpoint.index()];
            let explicit_hs = atom.explicit_hydrogens();
            if explicit_hs > 0 {
                let new_hs = (f64::from(explicit_hs) - order_difference) as i32;
                atom.set_explicit_hydrogens(new_hs.max(0) as u8);
            }
        }
    }
    if preserve_props {
        copied.replace_source_properties_from(old);
    }
    topology.bonds[bond.index()] = copied;
    if !keep_sgroups {
        crate::sgroup::remove_source_groups_referencing_bond(
            &mut topology.substance_groups,
            bond,
            &mut uint_reader,
        )?;
    }
    // BondId is unchanged: every occurrence in detached bookmarks already
    // refers to the replacement. Native stereo-group dangling pointers have
    // no stable-value equivalent; their modeled IDs remain unchanged.
    Ok(())
}

/// Canonical range-checked source edge lookup in original encounter order.
/// The callback is reached only after both atom range preconditions.
#[doc(hidden)]
pub fn source_bond_between_atoms<I: IntoIterator>(
    atom_count: usize,
    begin: AtomId,
    end: AtomId,
    neighbors: impl FnOnce() -> I,
) -> Result<Option<BondId>, TopologyEditError>
where
    I::Item: std::borrow::Borrow<crate::NeighborRef>,
{
    // RDKit❗✔️: const Bond *ROMol::getBondBetweenAtoms(unsigned int idx1,
    // RDKit❗✔️:                                        unsigned int idx2) const {
    // RDKit❗✔️:   URANGE_CHECK(idx1, getNumAtoms());
    // RDKit❗✔️:   URANGE_CHECK(idx2, getNumAtoms());
    // RDKit❗✔️:   const Bond *res = nullptr;
    // RDKit❗✔️:
    // RDKit❗✔️:   auto [edge, found] = boost::edge(boost::vertex(idx1, d_graph),
    // RDKit❗✔️:                                    boost::vertex(idx2, d_graph), d_graph);
    // RDKit❗✔️:   if (found) {
    // RDKit❗✔️:     res = d_graph[edge];
    // RDKit❗✔️:   }
    // RDKit❗✔️:   return res;
    // RDKit❗✔️: }
    for atom in [begin, end] {
        if atom.index() >= atom_count {
            return Err(TopologyEditError::AtomOutOfRange { atom, atom_count });
        }
    }
    Ok(neighbors()
        .into_iter()
        .find(|n| std::borrow::Borrow::<crate::NeighborRef>::borrow(n).atom_index == end.index())
        .map(|n| std::borrow::Borrow::<crate::NeighborRef>::borrow(&n).bond))
}

/// Canonical detached source bond-pointer overload. Owned input projects
/// takeOwnership=true; borrowed input is copied only after source preconditions.
#[doc(hidden)]
pub fn add_source_bond_value(
    atom_count: usize,
    bonds: &mut Vec<Bond>,
    neighbors: SourceBondNeighbors<'_>,
    bond: std::borrow::Cow<'_, Bond>,
) -> Result<usize, TopologyEditError> {
    // RDKit❗✔️: unsigned int ROMol::addBond(Bond *bond_pin, bool takeOwnership) {
    // RDKit❗✔️:   PRECONDITION(bond_pin, "null bond passed in");
    // RDKit❗✔️:   PRECONDITION(!takeOwnership || !bond_pin->hasOwningMol() ||
    // RDKit❗✔️:                    &bond_pin->getOwningMol() == this,
    // RDKit❗✔️:                "cannot take ownership of an bond which already has an owner");
    // RDKit❗✔️:   URANGE_CHECK(bond_pin->getBeginAtomIdx(), getNumAtoms());
    // RDKit❗✔️:   URANGE_CHECK(bond_pin->getEndAtomIdx(), getNumAtoms());
    // RDKit❗✔️:   PRECONDITION(bond_pin->getBeginAtomIdx() != bond_pin->getEndAtomIdx(),
    // RDKit❗✔️:                "attempt to add self-bond");
    // RDKit❗✔️:   PRECONDITION(!(boost::edge(bond_pin->getBeginAtomIdx(),
    // RDKit❗✔️:                              bond_pin->getEndAtomIdx(), d_graph)
    // RDKit❗✔️:                      .second),
    // RDKit❗✔️:                "bond already exists");
    // RDKit❗✔️:
    // RDKit❗✔️:   Bond *bond_p;
    // RDKit❗✔️:   if (!takeOwnership) {
    // RDKit❗✔️:     bond_p = bond_pin->copy();
    // RDKit❗✔️:   } else {
    // RDKit❗✔️:     bond_p = bond_pin;
    // RDKit❗✔️:   }
    // RDKit❗✔️:
    // RDKit❗✔️:   bond_p->setOwningMol(this);
    // RDKit❗✔️:   auto [which, ok] = boost::add_edge(bond_p->getBeginAtomIdx(),
    // RDKit❗✔️:                                      bond_p->getEndAtomIdx(), d_graph);
    // RDKit❗✔️:   CHECK_INVARIANT(ok, "bond could not be added");
    // RDKit❗✔️:   d_graph[which] = bond_p;
    // RDKit❗✔️:   bond_p->setIdx(numBonds);
    // RDKit❗✔️:   numBonds++;
    // RDKit❗✔️:   return numBonds;
    // RDKit❗✔️: }
    // Detached values have no null or live owner pointer. Both source ownership
    // choices are explicit; no atom/aromatic/cache/batch-mask effects are added.
    // Cost: same O(degree) edge check and amortized row/neighbor appends. Owned
    // trees move unchanged; borrowed trees clone exactly once after checks.
    let begin = bond.begin();
    let end = bond.end();
    for atom in [begin, end] {
        if atom.index() >= atom_count {
            return Err(TopologyEditError::AtomOutOfRange { atom, atom_count });
        }
    }
    if begin == end {
        return Err(TopologyEditError::InvalidResult(
            TopologyValidationError::SelfLoopBond {
                bond: BondId::new(bonds.len()),
                atom: begin,
            },
        ));
    }
    if neighbors.appended.len() != atom_count {
        return Err(TopologyEditError::InvalidSource(
            TopologyValidationError::AdjacencyMismatch,
        ));
    }
    if neighbors.original.is_some_and(|adj| {
        adj.neighbors_of(begin.index())
            .iter()
            .any(|n| n.atom_index == end.index())
    }) || neighbors.appended[begin.index()]
        .iter()
        .any(|n| n.atom_index == end.index())
    {
        return Err(TopologyEditError::DuplicateBond { begin, end });
    }
    let mut bond = bond.into_owned();
    let id = BondId::new(bonds.len());
    bond.set_id_for_construction(id);
    bonds.push(bond);
    neighbors.appended[begin.index()].push(crate::NeighborRef {
        atom_index: end.index(),
        bond: id,
    });
    neighbors.appended[end.index()].push(crate::NeighborRef {
        atom_index: begin.index(),
        bond: id,
    });
    Ok(bonds.len())
}

#[cfg(test)]
mod complete_source_add_bond_value_tests {
    use super::*;
    use crate::{
        AdjacencyList, BondDirection, BondOrder, BondQueryPredicate, BondSpec, PropertyValue,
        QueryNode,
    };
    use std::borrow::Cow;
    fn value() -> Bond {
        Bond::from_spec(
            BondId::new(17),
            BondSpec::new(AtomId::new(1), AtomId::new(0), BondOrder::Aromatic)
                .with_query(QueryNode::Predicate(BondQueryPredicate::Any))
                .with_direction(BondDirection::BeginWedge)
                .with_conjugated(true)
                .with_prop("source", PropertyValue::UInt(9))
                .unwrap(),
        )
    }
    #[test]
    fn owned_value_moves_query_without_clone_and_preserves_rich_bond_state() {
        let value = value();
        let pointer = value.query().unwrap() as *const _;
        let mut expected = value.clone();
        expected.set_id_for_construction(BondId::new(0));
        let mut bonds = vec![];
        let mut neighbors = vec![vec![]; 2];
        assert_eq!(
            add_source_bond_value(
                2,
                &mut bonds,
                SourceBondNeighbors {
                    original: None,
                    appended: &mut neighbors
                },
                Cow::Owned(value)
            )
            .unwrap(),
            1
        );
        assert_eq!(bonds, [expected]);
        assert_eq!(bonds[0].query().unwrap() as *const _, pointer);
        assert!(!bonds[0].is_aromatic());
        assert_eq!(neighbors[1][0].atom_index, 0);
        assert_eq!(neighbors[0][0].bond, BondId::new(0));
    }
    #[test]
    fn borrowed_value_copies_query_and_does_not_modify_source_id_or_properties() {
        let source = value();
        let before = source.clone();
        let mut bonds = vec![];
        let mut neighbors = vec![vec![]; 2];
        add_source_bond_value(
            2,
            &mut bonds,
            SourceBondNeighbors {
                original: None,
                appended: &mut neighbors,
            },
            Cow::Borrowed(&source),
        )
        .unwrap();
        assert_eq!(source, before);
        assert_eq!(bonds[0].query(), source.query());
        assert!(!std::ptr::eq(
            bonds[0].query().unwrap(),
            source.query().unwrap()
        ));
        assert_eq!(bonds[0].id(), BondId::new(0));
        assert_eq!(source.id(), BondId::new(17));
    }
    #[test]
    fn range_then_self_then_duplicate_preconditions_leave_all_rows_unchanged() {
        let mut bonds = vec![];
        let mut neighbors = vec![vec![]; 2];
        for (begin, end, kind) in [(99, 98, 0), (0, 98, 1), (0, 0, 2)] {
            let value = Bond::from_spec(
                BondId::new(9),
                BondSpec::new(AtomId::new(begin), AtomId::new(end), BondOrder::Single),
            );
            let error = add_source_bond_value(
                2,
                &mut bonds,
                SourceBondNeighbors {
                    original: None,
                    appended: &mut neighbors,
                },
                Cow::Owned(value),
            )
            .unwrap_err();
            match kind {
                0 => assert!(
                    matches!(error,TopologyEditError::AtomOutOfRange {atom,..} if atom==AtomId::new(99))
                ),
                1 => assert!(
                    matches!(error,TopologyEditError::AtomOutOfRange {atom,..} if atom==AtomId::new(98))
                ),
                _ => assert!(matches!(
                    error,
                    TopologyEditError::InvalidResult(TopologyValidationError::SelfLoopBond { .. })
                )),
            }
            assert!(bonds.is_empty());
            assert!(neighbors.iter().all(Vec::is_empty));
        }
        add_source_bond_value(
            2,
            &mut bonds,
            SourceBondNeighbors {
                original: None,
                appended: &mut neighbors,
            },
            Cow::Owned(value()),
        )
        .unwrap();
        let before = bonds.clone();
        let ns = neighbors.clone();
        assert!(matches!(
            add_source_bond_value(
                2,
                &mut bonds,
                SourceBondNeighbors {
                    original: None,
                    appended: &mut neighbors
                },
                Cow::Owned(value())
            ),
            Err(TopologyEditError::DuplicateBond { .. })
        ));
        assert_eq!(bonds, before);
        assert_eq!(neighbors, ns);
    }
    #[test]
    fn existing_source_adjacency_is_considered_before_append_and_owned_other_order_is_accepted() {
        let existing = vec![Bond::from_spec(
            BondId::new(0),
            BondSpec::new(AtomId::new(0), AtomId::new(1), BondOrder::Single),
        )];
        let adjacency = AdjacencyList::from_topology(3, &existing);
        let mut bonds = existing;
        let mut neighbors = vec![vec![]; 3];
        assert!(matches!(
            add_source_bond_value(
                3,
                &mut bonds,
                SourceBondNeighbors {
                    original: Some(&adjacency),
                    appended: &mut neighbors
                },
                Cow::Owned(value())
            ),
            Err(TopologyEditError::DuplicateBond { .. })
        ));
        let other = Bond::from_spec(
            BondId::new(9),
            BondSpec::new(AtomId::new(2), AtomId::new(0), BondOrder::Other),
        );
        assert_eq!(
            add_source_bond_value(
                3,
                &mut bonds,
                SourceBondNeighbors {
                    original: Some(&adjacency),
                    appended: &mut neighbors
                },
                Cow::Owned(other)
            )
            .unwrap(),
            2
        );
        assert_eq!(bonds[1].order(), BondOrder::Other);
        assert_eq!(bonds[1].id(), BondId::new(1));
    }
}

#[cfg(test)]
mod source_editor_commit_consumer_tests {
    use super::*;
    #[test]
    fn appended_edges_are_visible_to_endpoint_and_atom_removals_before_source_commit() {
        for atom_removal in [false, true] {
            let t = TopologyBlock::try_from_parts(
                (0..3)
                    .map(|i| Atom::from_spec(AtomId::new(i), AtomSpec::new(crate::Element::C)))
                    .collect(),
                vec![],
                vec![],
                vec![],
            )
            .unwrap();
            let mut e = t.into_batch_edit().unwrap();
            e.add_bond_by_order(AtomId::new(0), AtomId::new(2), crate::BondOrder::Single)
                .unwrap();
            if atom_removal {
                e.remove_atom(AtomId::new(2)).unwrap();
            } else {
                e.remove_bond_between_atoms(AtomId::new(0), AtomId::new(2))
                    .unwrap();
            }
            assert!(
                e.bond_between_atoms(AtomId::new(0), AtomId::new(2))
                    .unwrap()
                    .is_some()
            );
            let (t, m) = e
                .finish_source::<TopologyEditError>(
                    &mut Default::default(),
                    &mut Default::default(),
                    &mut |_| panic!("no SG properties"),
                    &mut |_| panic!("no ENDPTS"),
                )
                .unwrap();
            assert!(t.bonds.is_empty());
            assert_eq!(t.atoms.len(), if atom_removal { 2 } else { 3 });
            assert!(m.bonds.new_to_old.is_empty());
            assert!(m.bonds.old_to_new.is_empty());
            m.validate_for_counts(3, t.atoms.len(), 0, 0).unwrap();
        }
    }
}
