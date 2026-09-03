//! Coordinate and conformer value types shared by molecule algorithms.
//!
//! These values are detached working state.  The live `Molecule` owner and
//! topology/cache lifecycle remain in the runtime crate.

use std::collections::BTreeMap;

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

    #[allow(dead_code)]
    pub fn remapped_to_kept_atoms(&self, kept_old_indices: &[usize], id: usize) -> Self {
        let coords = kept_old_indices
            .iter()
            .filter_map(|old_idx| self.coords.get(*old_idx).copied())
            .collect();
        Self {
            id,
            coords,
            props: self.props.clone(),
        }
    }

    #[allow(dead_code)]
    pub fn push_coord(&mut self, coord: [f64; 2]) {
        self.coords.push(coord);
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
        self.id
    }

    #[must_use]
    pub fn coordinates(&self) -> &[[f64; 3]] {
        &self.coords
    }

    /// Mutable access to coordinates (pub(crate) for conformer transforms).
    pub fn coordinates_mut(&mut self) -> &mut [[f64; 3]] {
        &mut self.coords
    }

    #[must_use]
    pub const fn is_3d(&self) -> bool {
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
        self.id = id;
        self
    }

    #[allow(dead_code)]
    pub fn remapped_to_kept_atoms(&self, kept_old_indices: &[usize], id: usize) -> Self {
        let coords = kept_old_indices
            .iter()
            .filter_map(|old_idx| self.coords.get(*old_idx).copied())
            .collect();
        Self {
            id,
            coords,
            is_3d: self.is_3d,
            props: self.props.clone(),
        }
    }

    #[allow(dead_code)]
    pub fn push_coord(&mut self, coord: [f64; 3]) {
        self.coords.push(coord);
    }
}
#[derive(Debug, Clone, PartialEq, Default)]
pub struct ConformerStore {
    pub conformers_2d: Vec<Conformer2D>,
    pub conformers_3d: Vec<Conformer3D>,
    pub source_coordinate_dim: Option<CoordinateDimension>,
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
