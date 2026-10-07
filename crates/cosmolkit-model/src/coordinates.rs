//! Coordinate and conformer value types shared by molecule algorithms.
//!
//! These values are detached working state.  The live `Molecule` owner and
//! topology/cache lifecycle remain in the runtime crate.

use crate::PropertyText;
use std::collections::BTreeMap;

/// Locate the first non-finite value in row/column order for checked inputs.
/// Raw conformer storage and structural validation do not call this helper.
pub fn first_non_finite_coordinate<R: AsRef<[f64]>>(rows: &[R]) -> Option<(usize, usize, f64)> {
    for (row, values) in rows.iter().enumerate() {
        for (column, value) in values.as_ref().iter().copied().enumerate() {
            if !value.is_finite() {
                return Some((row, column, value));
            }
        }
    }
    None
}

#[derive(Debug, Clone, PartialEq, Eq, thiserror::Error)]
pub enum CoordinateValidationError {
    #[error("3D conformer id {id} does not exist")]
    Missing3DConformer { id: usize },
    #[error("3D conformer id overflow after {max_id}")]
    ConformerIdOverflow { max_id: usize },
    #[error("{dimension} conformer {conformer} has {rows} coordinate rows, expected {atom_count}")]
    RowCount {
        dimension: &'static str,
        conformer: usize,
        rows: usize,
        atom_count: usize,
    },
    #[error("mixed coordinate collections lack their source conformer order")]
    MissingSourceConformerOrder,
    #[error(
        "source conformer order contains {two_d} 2D and {three_d} 3D entries; expected {expected_two_d} and {expected_three_d}"
    )]
    SourceConformerOrder {
        two_d: usize,
        three_d: usize,
        expected_two_d: usize,
        expected_three_d: usize,
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
    props: BTreeMap<PropertyText, PropertyText>,
}

impl Conformer2D {
    /// Validate explicit checked input, preserving shape-before-number errors.
    pub fn validate_checked_for_atom_count(
        &self,
        atom_count: usize,
    ) -> Result<(), CoordinateValidationError> {
        self.validate_for_atom_count(atom_count)?;
        if let Some((atom, column, _)) = first_non_finite_coordinate(&self.coords) {
            return Err(CoordinateValidationError::NonFiniteCoordinate {
                dimension: "2D",
                conformer: self.id,
                atom,
                axis: ["x", "y"][column],
            });
        }
        Ok(())
    }

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
        // ROMol::addConformer checks the atom count, not numerical finiteness.
        // RDKit✔️✔️: PRECONDITION(conf->getNumAtoms() == this->getNumAtoms(),
        // RDKit✔️✔️:              "Number of atom mismatch");
        // Source Conformer::setAtomPos stores the supplied Point3D unchanged.
        // Non-finite values and their sign/payload are legitimate stored state;
        // any numerical algorithm's source-defined checks belong to that owner.
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
    pub fn props(&self) -> &BTreeMap<PropertyText, PropertyText> {
        &self.props
    }

    #[must_use]
    pub fn with_prop(
        mut self,
        key: impl Into<PropertyText>,
        value: impl Into<PropertyText>,
    ) -> Self {
        self.props.insert(key.into(), value.into());
        self
    }

    #[must_use]
    pub fn with_id(mut self, id: usize) -> Self {
        self.id = id;
        self
    }

    fn remapped_to_kept_atoms(&self, kept_old_indices: &[usize]) -> Self {
        let coords = kept_old_indices
            .iter()
            .map(|old_idx| self.coords[*old_idx])
            .collect();
        Self {
            id: self.id,
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
    props: BTreeMap<PropertyText, PropertyText>,
}

impl Conformer3D {
    /// Validate explicit checked input, preserving shape-before-number errors.
    pub fn validate_checked_for_atom_count(
        &self,
        atom_count: usize,
    ) -> Result<(), CoordinateValidationError> {
        self.validate_for_atom_count(atom_count)?;
        if let Some((atom, column, _)) = first_non_finite_coordinate(&self.coords) {
            return Err(CoordinateValidationError::NonFiniteCoordinate {
                dimension: "3D",
                conformer: self.id,
                atom,
                axis: ["x", "y", "z"][column],
            });
        }
        Ok(())
    }

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
        // ROMol::addConformer checks the atom count, not numerical finiteness.
        // RDKit✔️✔️: PRECONDITION(conf->getNumAtoms() == this->getNumAtoms(),
        // RDKit✔️✔️:              "Number of atom mismatch");
        // Source Conformer::setAtomPos stores the supplied Point3D unchanged.
        // Non-finite values and their sign/payload are legitimate stored state;
        // any numerical algorithm's source-defined checks belong to that owner.
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
    pub fn props(&self) -> &BTreeMap<PropertyText, PropertyText> {
        &self.props
    }

    #[must_use]
    pub fn with_prop(
        mut self,
        key: impl Into<PropertyText>,
        value: impl Into<PropertyText>,
    ) -> Self {
        self.props.insert(key.into(), value.into());
        self
    }

    #[must_use]
    pub fn with_id(mut self, id: usize) -> Self {
        // RDKit✔️✔️: inline void setId(unsigned int id) { d_id = id; }
        self.id = id;
        self
    }

    fn remapped_to_kept_atoms(&self, kept_old_indices: &[usize]) -> Self {
        let coords = kept_old_indices
            .iter()
            .map(|old_idx| self.coords[*old_idx])
            .collect();
        Self {
            id: self.id,
            coords,
            is_3d: self.is_3d,
            props: self.props.clone(),
        }
    }
}

#[derive(Debug, Clone, Default)]
pub struct CoordinateBlock {
    /// Zero or more 2D conformers.
    ///
    /// Coordinates are stored in the same atom-index order as `TopologyBlock`.
    /// Any operation changing atom indices must remap or drop this block through
    /// a topology report. Do not mutate this block directly from operation code.
    pub conformers_2d: Vec<Conformer2D>,
    pub conformers_3d: Vec<Conformer3D>,
    pub source_coordinate_dim: Option<CoordinateDimension>,
    /// The actual interleaving of source conformer appends. Occurrences index
    /// each dimension's collection; IDs and is_3d do not encode this fact.
    /// None is sufficient only when at most one dimension has stored rows.
    pub source_conformer_order: Option<Vec<CoordinateDimension>>,
}

impl PartialEq for CoordinateBlock {
    fn eq(&self, other: &Self) -> bool {
        if self.conformers_2d != other.conformers_2d
            || self.conformers_3d != other.conformers_3d
            || self.source_coordinate_dim != other.source_coordinate_dim
        {
            return false;
        }
        match (&self.source_conformer_order, &other.source_conformer_order) {
            (Some(left), Some(right)) => left == right,
            (None, None) => true,
            (Some(order), None) | (None, Some(order)) => {
                // A single typed collection already determines its complete
                // insertion order. Explicitly recording that same fact does
                // not change the detached coordinate value. Mixed unknown
                // order remains distinct from every explicit interleaving.
                if self.conformers_2d.is_empty() {
                    order.len() == self.conformers_3d.len()
                        && order.iter().all(|d| *d == CoordinateDimension::ThreeD)
                } else if self.conformers_3d.is_empty() {
                    order.len() == self.conformers_2d.len()
                        && order.iter().all(|d| *d == CoordinateDimension::TwoD)
                } else {
                    false
                }
            }
        }
    }
}

impl CoordinateBlock {
    pub fn validate_for_atom_count(
        &self,
        atom_count: usize,
    ) -> Result<(), CoordinateValidationError> {
        // RDKit✔️✔️: d_confs.push_back(nConf);
        // Source addConformer(assignId=false) does not deduplicate IDs.
        // Check every stored row while retaining vector order and all IDs.
        // Cost: one borrowed O(C) pass, without an extra ID set allocation.
        for conformer in &self.conformers_2d {
            conformer.validate_for_atom_count(atom_count)?;
        }
        for conformer in &self.conformers_3d {
            conformer.validate_for_atom_count(atom_count)?;
        }
        if let Some(order) = &self.source_conformer_order {
            let two_d = order
                .iter()
                .filter(|dimension| **dimension == CoordinateDimension::TwoD)
                .count();
            let three_d = order.len() - two_d;
            if two_d != self.conformers_2d.len() || three_d != self.conformers_3d.len() {
                return Err(CoordinateValidationError::SourceConformerOrder {
                    two_d,
                    three_d,
                    expected_two_d: self.conformers_2d.len(),
                    expected_three_d: self.conformers_3d.len(),
                });
            }
        }
        Ok(())
    }

    /// Remap conformer rows after a topology operation removes atoms.
    ///
    /// The runtime computes the authoritative topology mapping and passes only
    /// the retained old atom rows into this local value operation.
    /// Conformer IDs, collection order, properties and dimensionality are retained.
    pub fn remap_topology(&mut self, kept_old_indices: &[usize]) {
        // Pinned RDKit RWMol.cpp:901-913 (batchRemoveAtoms coordinate stage).
        // RDKit✔️❌: positions-only edit; existing conformer identity is unchanged.
        //   for (auto conf : d_confs) {
        //     RDGeom::POINT3D_VECT &positions = conf->getPositions();
        //     RDGeom::POINT3D_VECT newPositions;
        //     newPositions.reserve(getNumAtoms());
        //
        //     for (RDGeom::POINT3D_VECT::size_type i = 0; i < positions.size(); ++i) {
        //       if (oldIndices[i] != nullptr) {
        //         newPositions.push_back(positions[i]);
        //       }
        //     }
        //     CHECK_INVARIANT(newPositions.size() == getNumAtoms(), "Lost coordinates!");
        //     positions.swap(newPositions);
        //   }
        // Behavior review: project only atom rows, retaining each conformer's
        // ID/properties/is_3d and the block's source dimension. CK's authoritative
        // mapping additionally permits atom reordering; this is not a claim that
        // RDKit batch deletion reorders atoms. Mapping validation remains upstream.
        // Cost review: the existing implementation rebuilds both conformer Vecs,
        // allocates one coordinate Vec per conformer and clones its properties,
        // unlike the source's in-place conformer objects. This fix removes only
        // ID reassignment, not those existing allocation/copy costs.
        self.conformers_2d = self
            .conformers_2d
            .iter()
            .map(|conformer| conformer.remapped_to_kept_atoms(kept_old_indices))
            .collect();

        self.conformers_3d = self
            .conformers_3d
            .iter()
            .map(|conformer| conformer.remapped_to_kept_atoms(kept_old_indices))
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
    fn validate_preserves_duplicate_source_ids_per_dimension() {
        let block = CoordinateBlock {
            conformers_3d: vec![
                Conformer3D::new(4, vec![[0.0, 0.0, 0.0]], true),
                Conformer3D::new(4, vec![[1.0, 1.0, 1.0]], true),
            ],
            ..Default::default()
        };
        assert_eq!(block.validate_for_atom_count(1), Ok(()));
        assert_eq!(
            block
                .conformers_3d
                .iter()
                .map(Conformer3D::id)
                .collect::<Vec<_>>(),
            [4, 4]
        );
        assert_eq!(block.conformers_3d[1].coordinates(), &[[1.0, 1.0, 1.0]]);
    }
}

impl CoordinateBlock {
    /// Replace only the selected dimension-local ID, preserving row metadata.
    pub fn replace_3d_coordinates(
        &mut self,
        coords: Vec<[f64; 3]>,
        id: usize,
        atom_count: usize,
    ) -> Result<(), CoordinateValidationError> {
        let candidate = Conformer3D::new(id, coords, true);
        candidate.validate_for_atom_count(atom_count)?;
        let row = self
            .conformers_3d
            .iter_mut()
            .find(|row| row.id == id)
            .ok_or(CoordinateValidationError::Missing3DConformer { id })?;
        row.coords = candidate.coords;
        Ok(())
    }
    /// Clear the independent 3D dimension; 2D state is unchanged.
    pub fn clear_3d_conformers(&mut self) {
        // RDKit❗✔️: void clearConformers() { d_confs.clear(); }
        // Approved dimension-specific generation keeps every 2D row.
        self.conformers_3d.clear();
        if let Some(order) = &mut self.source_conformer_order {
            order.retain(|dimension| *dimension != CoordinateDimension::ThreeD);
        }
    }
    /// Install already-generated detached rows, without live commit authority.
    pub fn install_generated_3d(
        &mut self,
        clear_existing: bool,
        rows: Vec<Conformer3D>,
    ) -> Result<(), CoordinateValidationError> {
        if clear_existing {
            self.clear_3d_conformers();
        }
        for row in rows {
            self.record_source_conformer_append(CoordinateDimension::ThreeD)?;
            self.conformers_3d.push(row);
        }
        Ok(())
    }
    pub fn append_3d_conformer(
        &mut self,
        coords: Vec<[f64; 3]>,
        is_3d: bool,
        atom_count: usize,
        clear_existing: bool,
    ) -> Result<usize, CoordinateValidationError> {
        // RDKit❗✔️: unsigned int ROMol::addConformer(Conformer *conf, bool assignId) {
        // RDKit❗✔️:   PRECONDITION(conf, "bad conformer");
        // RDKit❗✔️:   PRECONDITION(conf->getNumAtoms() == this->getNumAtoms(),
        // RDKit❗✔️:                "Number of atom mismatch");
        // RDKit❗✔️:   if (assignId) {
        // RDKit❗✔️:     int maxId = -1;
        // RDKit❗✔️:     for (auto cptr : d_confs) {
        // RDKit❗✔️:       maxId = std::max((int)(cptr->getId()), maxId);
        // RDKit❗✔️:     }
        // RDKit❗✔️:     maxId++;
        // RDKit❗✔️:     conf->setId((unsigned int)maxId);
        // RDKit❗✔️:   }
        // RDKit❗✔️:   conf->setOwningMol(this);
        // RDKit❗✔️:   CONFORMER_SPTR nConf(conf);
        // RDKit❗✔️:   d_confs.push_back(nConf);
        // RDKit❗✔️:   return conf->getId();
        // RDKit❗✔️: }
        // Approved typed collections scope this source scan to 3D. O(C) borrowed scan;
        // validate first, then move one row. Wider identifier overflow is typed.
        let id = if clear_existing {
            0
        } else {
            match self.conformers_3d.iter().map(Conformer3D::id).max() {
                None => 0,
                Some(max_id) => max_id
                    .checked_add(1)
                    .ok_or(CoordinateValidationError::ConformerIdOverflow { max_id })?,
            }
        };
        let candidate = Conformer3D::new(id, coords, is_3d);
        candidate.validate_for_atom_count(atom_count)?;
        if clear_existing {
            self.clear_3d_conformers();
        }
        self.record_source_conformer_append(CoordinateDimension::ThreeD)?;
        self.conformers_3d.push(candidate);
        Ok(id)
    }
}

/// A borrowed coordinate set in its actual source insertion order.
#[derive(Debug, Clone, Copy)]
pub enum CoordinateSourceConformer<'a> {
    TwoD(&'a Conformer2D),
    ThreeD(&'a Conformer3D),
}

impl CoordinateBlock {
    /// Return the source's first conformer, without consulting IDs, numerical
    /// geometry, dimensional flags, or import-provenance hints.
    pub fn first_source_conformer(
        &self,
    ) -> Result<Option<CoordinateSourceConformer<'_>>, CoordinateValidationError> {
        first_source_conformer_from_parts(
            &self.conformers_2d,
            &self.conformers_3d,
            self.source_conformer_order.as_deref(),
        )
    }

    /// Record an actual append, before its row is pushed into typed storage.
    pub fn record_source_conformer_append(
        &mut self,
        dimension: CoordinateDimension,
    ) -> Result<(), CoordinateValidationError> {
        // RDKit❗✔️: d_confs.push_back(nConf);
        // Once materialized, a dimension occurrence appends in amortized O(1).
        // A previously dimension-only detached collection has an unambiguous
        // existing order; materialization costs O(C) only on the first append.
        record_source_append(
            &mut self.source_conformer_order,
            self.conformers_2d.len(),
            self.conformers_3d.len(),
            dimension,
        )?;
        Ok(())
    }

    /// Remove one typed dimension while retaining every other occurrence's
    /// relative source order. Replacing rows must record subsequent appends.
    pub fn clear_2d_conformers(&mut self) {
        self.conformers_2d.clear();
        if let Some(order) = &mut self.source_conformer_order {
            order.retain(|dimension| *dimension != CoordinateDimension::TwoD);
        }
    }
}

pub(crate) fn record_source_append(
    order: &mut Option<Vec<CoordinateDimension>>,
    two_d: usize,
    three_d: usize,
    dimension: CoordinateDimension,
) -> Result<(), CoordinateValidationError> {
    // RDKit❗✔️: d_confs.push_back(nConf);
    // Appending does not inspect an existing conformer's ID or geometry.
    // A detached mixed prefix without its source order cannot recover that
    // fact by appending. Preserve its explicit unknown state instead of
    // guessing an order or rejecting a dimension-local append. A source-first
    // read still returns MissingSourceConformerOrder for that prefix.
    // Known order appends in amortized O(1); a single-dimension prefix is
    // materialized once in O(C), without copying any coordinate rows.
    if order.is_none() {
        let (existing_dimension, count) = if two_d == 0 {
            (CoordinateDimension::ThreeD, three_d)
        } else if three_d == 0 {
            (CoordinateDimension::TwoD, two_d)
        } else {
            return Ok(());
        };
        *order = Some(vec![existing_dimension; count]);
    }
    order.as_mut().unwrap().push(dimension);
    Ok(())
}

pub(crate) fn first_source_conformer_from_parts<'a>(
    conformers_2d: &'a [Conformer2D],
    conformers_3d: &'a [Conformer3D],
    source_conformer_order: Option<&[CoordinateDimension]>,
) -> Result<Option<CoordinateSourceConformer<'a>>, CoordinateValidationError> {
    // RDKit❗✔️: const Conformer &ROMol::getConformer(int id) const {
    // RDKit❗✔️:   if (d_confs.size() == 0) {
    // RDKit❗✔️:     throw ConformerException("No conformations available on the molecule");
    // RDKit❗✔️:   }
    // RDKit❗✔️:   if (id < 0) {
    // RDKit❗✔️:     return *(d_confs.front());
    // RDKit❗✔️:   }
    // The detached optional empty result is translated by each caller's
    // source-defined empty branch. Nonempty lookup is O(1), like front().
    let first = match source_conformer_order {
        Some(order) => order.first().copied(),
        None if conformers_2d.is_empty() => {
            conformers_3d.first().map(|_| CoordinateDimension::ThreeD)
        }
        None if conformers_3d.is_empty() => {
            conformers_2d.first().map(|_| CoordinateDimension::TwoD)
        }
        None => return Err(CoordinateValidationError::MissingSourceConformerOrder),
    };
    Ok(match first {
        Some(CoordinateDimension::TwoD) => {
            conformers_2d.first().map(CoordinateSourceConformer::TwoD)
        }
        Some(CoordinateDimension::ThreeD) => {
            conformers_3d.first().map(CoordinateSourceConformer::ThreeD)
        }
        None => None,
    })
}
#[cfg(test)]
mod source_order_tests {
    use super::*;

    #[test]
    fn append_to_unknown_mixed_prefix_preserves_unknown_without_guessing() {
        let mut block = CoordinateBlock {
            conformers_2d: vec![Conformer2D::new(9, vec![[1.0, -0.0]])],
            conformers_3d: vec![Conformer3D::new(9, vec![[2.0, 3.0, 4.0]], false)],
            ..Default::default()
        };
        block
            .record_source_conformer_append(CoordinateDimension::ThreeD)
            .unwrap();
        block
            .conformers_3d
            .push(Conformer3D::new(0, vec![[5.0, 6.0, 7.0]], true));
        block.validate_for_atom_count(1).unwrap();
        assert_eq!(block.source_conformer_order, None);
        assert!(matches!(
            block.first_source_conformer(),
            Err(CoordinateValidationError::MissingSourceConformerOrder)
        ));
        assert_eq!(
            block.conformers_2d[0].coordinates()[0][1].to_bits(),
            (-0.0_f64).to_bits()
        );
        assert_eq!(block.conformers_3d[0].id(), 9);
        assert!(!block.conformers_3d[0].is_3d());
    }

    #[test]
    fn source_order_keeps_actual_first_across_dimensions_and_duplicate_ids() {
        let mut block = CoordinateBlock::default();
        block
            .record_source_conformer_append(CoordinateDimension::ThreeD)
            .unwrap();
        block
            .conformers_3d
            .push(Conformer3D::new(9, vec![[3.0, 4.0, 5.0]], false));
        block
            .record_source_conformer_append(CoordinateDimension::TwoD)
            .unwrap();
        block
            .conformers_2d
            .push(Conformer2D::new(9, vec![[1.0, 2.0]]));
        block
            .record_source_conformer_append(CoordinateDimension::ThreeD)
            .unwrap();
        block
            .conformers_3d
            .push(Conformer3D::new(0, vec![[6.0, 7.0, 8.0]], true));
        block.validate_for_atom_count(1).unwrap();
        match block.first_source_conformer().unwrap().unwrap() {
            CoordinateSourceConformer::ThreeD(conformer) => {
                assert_eq!(conformer.id(), 9);
                assert!(!conformer.is_3d());
                assert_eq!(conformer.coordinates(), &[[3.0, 4.0, 5.0]]);
            }
            _ => panic!("the first actual source append is retained"),
        }
        block.clear_3d_conformers();
        match block.first_source_conformer().unwrap().unwrap() {
            CoordinateSourceConformer::TwoD(conformer) => assert_eq!(conformer.id(), 9),
            _ => panic!("dimension-local clear retains remaining source order"),
        }
    }

    #[test]
    fn source_order_extent_validation_and_row_remap_preserve_facts() {
        let mut block = CoordinateBlock {
            conformers_2d: vec![Conformer2D::new(8, vec![[0.0, 1.0], [2.0, 3.0]])],
            conformers_3d: vec![Conformer3D::new(
                8,
                vec![[4.0, 5.0, 6.0], [7.0, 8.0, 9.0]],
                true,
            )],
            source_conformer_order: Some(vec![
                CoordinateDimension::ThreeD,
                CoordinateDimension::TwoD,
            ]),
            ..Default::default()
        };
        block.validate_for_atom_count(2).unwrap();
        block.remap_topology(&[1]);
        block.validate_for_atom_count(1).unwrap();
        assert_eq!(
            block.source_conformer_order,
            Some(vec![CoordinateDimension::ThreeD, CoordinateDimension::TwoD])
        );
        assert_eq!(block.conformers_3d[0].coordinates(), &[[7.0, 8.0, 9.0]]);
        block.source_conformer_order = Some(vec![CoordinateDimension::TwoD]);
        assert!(matches!(
            block.validate_for_atom_count(1),
            Err(CoordinateValidationError::SourceConformerOrder { .. })
        ));
    }
}
