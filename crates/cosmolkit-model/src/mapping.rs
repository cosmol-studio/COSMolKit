//! Detached atom and bond row mappings.
//!
//! This module validates local old/new index-space relationships. Applying a
//! mapping to live molecule state and authorizing an operation commit remain
//! responsibilities of the parent runtime.

use crate::{AtomId, BondId};

#[derive(Debug, Clone, PartialEq, Eq, thiserror::Error)]
pub enum MappingValidationError {
    #[error("{entity} {direction} mapping has {actual} rows, expected {expected}")]
    Length {
        entity: &'static str,
        direction: &'static str,
        actual: usize,
        expected: usize,
    },
    #[error(
        "{entity} {direction} mapping row {row} refers to {mapped}, outside {target_count} rows"
    )]
    OutOfRange {
        entity: &'static str,
        direction: &'static str,
        row: usize,
        mapped: usize,
        target_count: usize,
    },
    #[error("{entity} mappings disagree for {direction} row {row} and mapped row {mapped}")]
    InverseMismatch {
        entity: &'static str,
        direction: &'static str,
        row: usize,
        mapped: usize,
    },
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

    /// Return old atom rows in new-row order, omitting appended rows.
    ///
    /// The complete mapping must be validated against its authoritative old
    /// and new counts before this projection is passed to a remap helper.
    #[must_use]
    pub fn retained_atom_indices(&self) -> Vec<usize> {
        self.atoms
            .new_to_old
            .iter()
            .filter_map(|id| id.map(AtomId::index))
            .collect()
    }

    #[must_use]
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

    #[must_use]
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

    /// Validate both directions against the authoritative old and new table
    /// sizes before any dependent-state projection or runtime commit.
    pub fn validate_for_counts(
        &self,
        old_atom_count: usize,
        new_atom_count: usize,
        old_bond_count: usize,
        new_bond_count: usize,
    ) -> Result<(), MappingValidationError> {
        validate_atom_mapping(&self.atoms, old_atom_count, new_atom_count)?;
        validate_bond_mapping(&self.bonds, old_bond_count, new_bond_count)
    }
}

fn validate_atom_mapping(
    mapping: &AtomMapping,
    old_count: usize,
    new_count: usize,
) -> Result<(), MappingValidationError> {
    if mapping.old_to_new.len() != old_count {
        return Err(MappingValidationError::Length {
            entity: "atom",
            direction: "old-to-new",
            actual: mapping.old_to_new.len(),
            expected: old_count,
        });
    }
    if mapping.new_to_old.len() != new_count {
        return Err(MappingValidationError::Length {
            entity: "atom",
            direction: "new-to-old",
            actual: mapping.new_to_old.len(),
            expected: new_count,
        });
    }
    for (old, new) in mapping.old_to_new.iter().enumerate() {
        let Some(new) = new else {
            continue;
        };
        if new.index() >= new_count {
            return Err(MappingValidationError::OutOfRange {
                entity: "atom",
                direction: "old-to-new",
                row: old,
                mapped: new.index(),
                target_count: new_count,
            });
        }
        if mapping.new_to_old[new.index()] != Some(AtomId::new(old)) {
            return Err(MappingValidationError::InverseMismatch {
                entity: "atom",
                direction: "old-to-new",
                row: old,
                mapped: new.index(),
            });
        }
    }
    for (new, old) in mapping.new_to_old.iter().enumerate() {
        let Some(old) = old else {
            continue;
        };
        if old.index() >= old_count {
            return Err(MappingValidationError::OutOfRange {
                entity: "atom",
                direction: "new-to-old",
                row: new,
                mapped: old.index(),
                target_count: old_count,
            });
        }
        if mapping.old_to_new[old.index()] != Some(AtomId::new(new)) {
            return Err(MappingValidationError::InverseMismatch {
                entity: "atom",
                direction: "new-to-old",
                row: new,
                mapped: old.index(),
            });
        }
    }
    Ok(())
}

fn validate_bond_mapping(
    mapping: &BondMapping,
    old_count: usize,
    new_count: usize,
) -> Result<(), MappingValidationError> {
    if mapping.old_to_new.len() != old_count {
        return Err(MappingValidationError::Length {
            entity: "bond",
            direction: "old-to-new",
            actual: mapping.old_to_new.len(),
            expected: old_count,
        });
    }
    if mapping.new_to_old.len() != new_count {
        return Err(MappingValidationError::Length {
            entity: "bond",
            direction: "new-to-old",
            actual: mapping.new_to_old.len(),
            expected: new_count,
        });
    }
    for (old, new) in mapping.old_to_new.iter().enumerate() {
        let Some(new) = new else {
            continue;
        };
        if new.index() >= new_count {
            return Err(MappingValidationError::OutOfRange {
                entity: "bond",
                direction: "old-to-new",
                row: old,
                mapped: new.index(),
                target_count: new_count,
            });
        }
        if mapping.new_to_old[new.index()] != Some(BondId::new(old)) {
            return Err(MappingValidationError::InverseMismatch {
                entity: "bond",
                direction: "old-to-new",
                row: old,
                mapped: new.index(),
            });
        }
    }
    for (new, old) in mapping.new_to_old.iter().enumerate() {
        let Some(old) = old else {
            continue;
        };
        if old.index() >= old_count {
            return Err(MappingValidationError::OutOfRange {
                entity: "bond",
                direction: "new-to-old",
                row: new,
                mapped: old.index(),
                target_count: old_count,
            });
        }
        if mapping.old_to_new[old.index()] != Some(BondId::new(new)) {
            return Err(MappingValidationError::InverseMismatch {
                entity: "bond",
                direction: "new-to-old",
                row: new,
                mapped: old.index(),
            });
        }
    }
    Ok(())
}
