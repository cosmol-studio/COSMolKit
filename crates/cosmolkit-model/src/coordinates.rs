//! Coordinate and conformer value types shared by molecule algorithms.
//!
//! These values are detached working state.  The live `Molecule` owner and
//! topology/cache lifecycle remain in the runtime crate.

use std::collections::BTreeMap;

#[derive(Debug, Clone, PartialEq, Eq, thiserror::Error)]
pub enum CoordinateValidationError {
    #[error("{dimension} conformer {conformer} has {rows} coordinate rows, expected {atom_count}")]
    RowCount {
        dimension: &'static str,
        conformer: usize,
        rows: usize,
        atom_count: usize,
    },
    #[error("duplicate {dimension} conformer id {id}")]
    DuplicateConformerId { dimension: &'static str, id: usize },
    #[error("{dimension} conformer {conformer} atom row {atom} has a non-finite {axis} coordinate")]
    NonFiniteCoordinate {
        dimension: &'static str,
        conformer: usize,
        atom: usize,
        axis: &'static str,
    },
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CoordinateDimension {
    TwoD,
    ThreeD,
}

#[derive(Debug, Clone, PartialEq)]
pub struct Conformer2D {
    id: usize,
    coords: Vec<[f64; 2]>,
    props: BTreeMap<String, String>,
}

impl Conformer2D {
    pub fn validate_for_atom_count(
        &self,
        atom_count: usize,
    ) -> Result<(), CoordinateValidationError> {
        if self.coords.len() != atom_count {
            return Err(CoordinateValidationError::RowCount {
                dimension: "2D",
                conformer: self.id,
                rows: self.coords.len(),
                atom_count,
            });
        }
        for (atom, coord) in self.coords.iter().enumerate() {
            for (axis, value) in [("x", coord[0]), ("y", coord[1])] {
                if !value.is_finite() {
                    return Err(CoordinateValidationError::NonFiniteCoordinate {
                        dimension: "2D",
                        conformer: self.id,
                        atom,
                        axis,
                    });
                }
            }
        }
        Ok(())
    }

    #[must_use]
    pub fn new(id: usize, coords: Vec<[f64; 2]>) -> Self {
        Self {
            id,
            coords,
            props: BTreeMap::new(),
        }
    }

    #[must_use]
    pub const fn id(&self) -> usize {
        self.id
    }

    #[must_use]
    pub fn coordinates(&self) -> &[[f64; 2]] {
        &self.coords
    }

    pub fn coordinates_mut(&mut self) -> &mut [[f64; 2]] {
        &mut self.coords
    }

    #[must_use]
    pub fn props(&self) -> &BTreeMap<String, String> {
        &self.props
    }

    #[must_use]
    pub fn with_prop(mut self, key: impl Into<String>, value: impl Into<String>) -> Self {
        self.props.insert(key.into(), value.into());
        self
    }

    #[must_use]
    pub fn with_id(mut self, id: usize) -> Self {
        self.id = id;
        self
    }

    fn remapped_to_kept_atoms(&self, kept_old_indices: &[usize], id: usize) -> Self {
        let coords = kept_old_indices
            .iter()
            .map(|old_idx| self.coords[*old_idx])
            .collect();
        Self {
            id,
            coords,
            props: self.props.clone(),
        }
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct Conformer3D {
    id: usize,
    coords: Vec<[f64; 3]>,
    is_3d: bool,
    props: BTreeMap<String, String>,
}

impl Conformer3D {
    pub fn validate_for_atom_count(
        &self,
        atom_count: usize,
    ) -> Result<(), CoordinateValidationError> {
        if self.coords.len() != atom_count {
            return Err(CoordinateValidationError::RowCount {
                dimension: "3D",
                conformer: self.id,
                rows: self.coords.len(),
                atom_count,
            });
        }
        for (atom, coord) in self.coords.iter().enumerate() {
            for (axis, value) in [("x", coord[0]), ("y", coord[1]), ("z", coord[2])] {
                if !value.is_finite() {
                    return Err(CoordinateValidationError::NonFiniteCoordinate {
                        dimension: "3D",
                        conformer: self.id,
                        atom,
                        axis,
                    });
                }
            }
        }
        Ok(())
    }

    #[must_use]
    pub fn new(id: usize, coords: Vec<[f64; 3]>, is_3d: bool) -> Self {
        Self {
            id,
            coords,
            is_3d,
            props: BTreeMap::new(),
        }
    }

    #[must_use]
    pub const fn id(&self) -> usize {
        // RDKit✔️✔️: inline unsigned int getId() const { return d_id; }
        self.id
    }

    #[must_use]
    pub fn coordinates(&self) -> &[[f64; 3]] {
        // RDKit✔️✔️: const RDGeom::POINT3D_VECT &Conformer::getPositions() const {
        // RDKit✔️✔️:   if (dp_mol) {
        // RDKit✔️✔️:     PRECONDITION(dp_mol->getNumAtoms() == d_positions.size(), "");
        // RDKit✔️✔️:   }
        // RDKit✔️✔️:   return d_positions;
        // RDKit✔️✔️: }
        // Detached model conformers have no owning-molecule pointer, so the
        // source ownership precondition is unreachable in the modeled state.
        &self.coords
    }

    pub fn coordinates_mut(&mut self) -> &mut [[f64; 3]] {
        // RDKit❗✔️: RDGeom::POINT3D_VECT &Conformer::getPositions() { return d_positions; }
        // COSMolKit intentionally exposes a fixed-length mutable slice here;
        // element mutation matches the source, while vector resizing is not
        // reproduced by this detached value accessor.
        &mut self.coords
    }

    #[must_use]
    pub const fn is_3d(&self) -> bool {
        // RDKit✔️✔️: inline bool is3D() const { return df_is3D; }
        self.is_3d
    }

    #[must_use]
    pub fn props(&self) -> &BTreeMap<String, String> {
        &self.props
    }

    #[must_use]
    pub fn with_prop(mut self, key: impl Into<String>, value: impl Into<String>) -> Self {
        self.props.insert(key.into(), value.into());
        self
    }

    #[must_use]
    pub fn with_id(mut self, id: usize) -> Self {
        // RDKit✔️✔️: inline void setId(unsigned int id) { d_id = id; }
        self.id = id;
        self
    }

    fn remapped_to_kept_atoms(&self, kept_old_indices: &[usize], id: usize) -> Self {
        let coords = kept_old_indices
            .iter()
            .map(|old_idx| self.coords[*old_idx])
            .collect();
        Self {
            id,
            coords,
            is_3d: self.is_3d,
            props: self.props.clone(),
        }
    }
}

#[derive(Debug, Clone, PartialEq, Default)]
pub struct CoordinateBlock {
    /// Zero or more 2D conformers.
    ///
    /// Coordinates are stored in the same atom-index order as `TopologyBlock`.
    /// Any operation changing atom indices must remap or drop this block through
    /// a topology report. Do not mutate this block directly from operation code.
    pub conformers_2d: Vec<Conformer2D>,
    pub conformers_3d: Vec<Conformer3D>,
    pub source_coordinate_dim: Option<CoordinateDimension>,
}

impl CoordinateBlock {
    pub fn validate_for_atom_count(
        &self,
        atom_count: usize,
    ) -> Result<(), CoordinateValidationError> {
        let mut ids = std::collections::BTreeSet::new();
        for conformer in &self.conformers_2d {
            if !ids.insert(conformer.id()) {
                return Err(CoordinateValidationError::DuplicateConformerId {
                    dimension: "2D",
                    id: conformer.id(),
                });
            }
            conformer.validate_for_atom_count(atom_count)?;
        }
        ids.clear();
        for conformer in &self.conformers_3d {
            if !ids.insert(conformer.id()) {
                return Err(CoordinateValidationError::DuplicateConformerId {
                    dimension: "3D",
                    id: conformer.id(),
                });
            }
            conformer.validate_for_atom_count(atom_count)?;
        }
        Ok(())
    }

    /// Remap conformer rows after a topology operation removes atoms.
    ///
    /// The runtime computes the authoritative topology mapping and passes only
    /// the retained old atom rows into this local value operation.
    pub fn remap_topology(&mut self, kept_old_indices: &[usize]) {
        self.conformers_2d = self
            .conformers_2d
            .iter()
            .enumerate()
            .map(|(id, conformer)| conformer.remapped_to_kept_atoms(kept_old_indices, id))
            .collect();

        self.conformers_3d = self
            .conformers_3d
            .iter()
            .enumerate()
            .map(|(id, conformer)| conformer.remapped_to_kept_atoms(kept_old_indices, id))
            .collect();
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn validate_rejects_coordinate_row_mismatch() {
        let block = CoordinateBlock {
            conformers_2d: vec![Conformer2D::new(0, vec![[0.0, 0.0]])],
            ..Default::default()
        };
        assert!(matches!(
            block.validate_for_atom_count(2),
            Err(CoordinateValidationError::RowCount {
                dimension: "2D",
                ..
            })
        ));
    }

    #[test]
    fn validate_rejects_duplicate_ids_per_dimension() {
        let block = CoordinateBlock {
            conformers_3d: vec![
                Conformer3D::new(4, vec![[0.0, 0.0, 0.0]], true),
                Conformer3D::new(4, vec![[1.0, 1.0, 1.0]], true),
            ],
            ..Default::default()
        };
        assert!(matches!(
            block.validate_for_atom_count(1),
            Err(CoordinateValidationError::DuplicateConformerId {
                dimension: "3D",
                id: 4
            })
        ));
    }
}
