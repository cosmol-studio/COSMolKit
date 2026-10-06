//! Project-native explicit coordinate installation from pinned COSMolKit d892.
//! Numeric matrix interpretation and detached transitions have this sole owner.
use cosmolkit_model::{
    Conformer2D, Conformer3D, CoordinateBlock, CoordinateDimension, CoordinateValidationError,
};

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum CoordinateZPolicy {
    #[default]
    Ignore,
    RequireZero,
    Error,
}
impl CoordinateZPolicy {
    pub fn from_name(value: &str) -> Result<Self, CoordinateInputError> {
        // Pinned COSMolKit d892 original numeric ingress, verbatim source:
        // fn parse_z_policy(value: &str) -> PyResult<&'static str> {
        //     match value.to_ascii_lowercase().as_str() {
        //         "ignore" => Ok("ignore"),
        //         "require_zero" => Ok("require_zero"),
        //         "error" => Ok("error"),
        //         _ => Err(PyValueError::new_err(format!(
        //             "unsupported z_policy '{value}', expected one of: ignore, require_zero, error"
        //         ))),
        //     }
        // }
        // COSMolKit✔️✔️: Original validation order and linear numeric conversion.

        // Pinned d892 parse_z_policy uses ASCII case folding only.
        match value.to_ascii_lowercase().as_str() {
            "ignore" => Ok(Self::Ignore),
            "require_zero" => Ok(Self::RequireZero),
            "error" => Ok(Self::Error),
            _ => Err(CoordinateInputError::UnknownZPolicy {
                value: value.to_owned(),
            }),
        }
    }
}
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct Coordinate2DInputParams {
    pub z_policy: CoordinateZPolicy,
}
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Coordinate3DInputParams {
    pub is_3d: bool,
}
impl Default for Coordinate3DInputParams {
    fn default() -> Self {
        Self { is_3d: true }
    }
}
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct Replace3DCoordinatesParams {
    /// Actual ThreeD conformer ID; the default selects ID0 and never falls back.
    pub conformer_id: usize,
}

#[derive(Clone, Debug, PartialEq, thiserror::Error)]
pub enum CoordinateInputError {
    #[error("unsupported z_policy '{value}', expected one of: ignore, require_zero, error")]
    UnknownZPolicy { value: String },
    #[error("{dimension} coordinates row count mismatch: expected {expected}, got {actual}")]
    RowCount {
        dimension: &'static str,
        expected: usize,
        actual: usize,
    },
    #[error(
        "{dimension} coordinates must have shape (num_atoms, {expected}); row {row} has {columns} columns"
    )]
    Shape {
        dimension: &'static str,
        row: usize,
        columns: usize,
        expected: &'static str,
    },
    #[error("{dimension} coordinates contains a non-finite value at row {row}, column {column}")]
    NonFinite {
        dimension: &'static str,
        row: usize,
        column: usize,
        value: f64,
    },
    #[error("2D coordinates require zero z values but row {row} has z={z}")]
    NonZeroZ { row: usize, z: f64 },
    #[error("2D coordinates received 3 columns while z_policy='error'")]
    ThreeColumnsForbidden,
    #[error("no 3D conformer with id {conformer_id}; {count} conformers stored")]
    ConformerNotFound { conformer_id: usize, count: usize },
    #[error("3D conformer identifier overflow after {max_id}")]
    ConformerIdOverflow { max_id: usize },
    #[error(transparent)]
    InvalidCoordinates(#[from] CoordinateValidationError),
}

fn validate_matrix(
    atom_count: usize,
    rows: &[Vec<f64>],
    dimension: &'static str,
    allowed: &[usize],
    expected: &'static str,
) -> Result<(), CoordinateInputError> {
    // Pinned COSMolKit d892 original numeric ingress, verbatim source:
    // fn extract_coordinate_matrix(
    //     coords: &Bound<'_, PyAny>,
    //     expected_rows: usize,
    //     expected_columns: &[usize],
    //     label: &str,
    // ) -> PyResult<Vec<Vec<f64>>> {
    //     if let Ok(array) = coords.cast::<PyUntypedArray>()
    //         && let [actual_rows, _] = array.shape()
    //         && *actual_rows != expected_rows
    //     {
    //         return Err(PyValueError::new_err(format!(
    //             "{label} row count mismatch: expected {expected_rows}, got {actual_rows}"
    //         )));
    //     }
    //     let array_like = coords
    //         .extract::<PyArrayLike<'_, f64, Ix2, AllowTypeChange>>()
    //         .map_err(|err| {
    //             PyTypeError::new_err(format!("{label} must be a 2D numeric array: {err}"))
    //         })?;
    //     let array = array_like.as_array();
    //     let shape = array.shape();
    //     if shape[0] != expected_rows {
    //         return Err(PyValueError::new_err(format!(
    //             "{label} row count mismatch: expected {expected_rows}, got {}",
    //             shape[0]
    //         )));
    //     }
    //     if !expected_columns.contains(&shape[1]) {
    //         let columns = expected_columns
    //             .iter()
    //             .map(ToString::to_string)
    //             .collect::<Vec<_>>()
    //             .join(" or ");
    //         return Err(PyValueError::new_err(format!(
    //             "{label} must have shape (num_atoms, {columns}); got ({}, {})",
    //             shape[0], shape[1]
    //         )));
    //     }
    //     let mut out = Vec::with_capacity(shape[0]);
    //     for (row_idx, row) in array.outer_iter().enumerate() {
    //         let mut values = Vec::with_capacity(shape[1]);
    //         for (col_idx, value) in row.iter().enumerate() {
    //             if !value.is_finite() {
    //                 return Err(PyValueError::new_err(format!(
    //                     "{label} contains a non-finite value at row {row_idx}, column {col_idx}"
    //                 )));
    //             }
    //             values.push(*value);
    //         }
    //         out.push(values);
    //     }
    //     Ok(out)
    // }
    // COSMolKit✔️✔️: Original validation order and linear numeric conversion.

    // Pinned d892 Python extract_coordinate_matrix: row count precedes shape;
    // every value, including subsequently discarded z, must be finite.
    // Behavior: original rectangular Python matrices keep their validation order.
    // Rust also validates each row of its explicit row representation.
    // Cost: one O(rows*columns) borrowed pass, without allocating/cloning rows.
    if rows.len() != atom_count {
        return Err(CoordinateInputError::RowCount {
            dimension,
            expected: atom_count,
            actual: rows.len(),
        });
    }
    for (row, values) in rows.iter().enumerate() {
        if !allowed.contains(&values.len()) {
            return Err(CoordinateInputError::Shape {
                dimension,
                row,
                columns: values.len(),
                expected,
            });
        }
    }
    if let Some((row, column, value)) = cosmolkit_model::first_non_finite_coordinate(rows) {
        return Err(CoordinateInputError::NonFinite {
            dimension,
            row,
            column,
            value,
        });
    }
    Ok(())
}

/// Interpret the complete original two/three-column input before any edit.
pub fn coordinates_2d_from_input(
    atom_count: usize,
    rows: Vec<Vec<f64>>,
    params: &Coordinate2DInputParams,
) -> Result<Vec<[f64; 2]>, CoordinateInputError> {
    // Pinned COSMolKit d892 original numeric ingress, verbatim source:
    // fn extract_2d_coordinates(
    //     coords: &Bound<'_, PyAny>,
    //     expected_rows: usize,
    //     z_policy: &str,
    // ) -> PyResult<Vec<[f64; 2]>> {
    //     let z_policy = parse_z_policy(z_policy)?;
    //     let rows = extract_coordinate_matrix(coords, expected_rows, &[2, 3], "2D coordinates")?;
    //     let mut out = Vec::with_capacity(rows.len());
    //     for (row_idx, row) in rows.into_iter().enumerate() {
    //         if row.len() == 3 {
    //             match z_policy {
    //                 "ignore" => {}
    //                 "require_zero" if row[2].abs() <= 1.0e-12 => {}
    //                 "require_zero" => {
    //                     return Err(PyValueError::new_err(format!(
    //                         "2D coordinates require zero z values but row {row_idx} has z={}",
    //                         row[2]
    //                     )));
    //                 }
    //                 "error" => {
    //                     return Err(PyValueError::new_err(
    //                         "2D coordinates received 3 columns while z_policy='error'",
    //                     ));
    //                 }
    //                 _ => unreachable!("z_policy is validated above"),
    //             }
    //         }
    //         out.push([row[0], row[1]]);
    //     }
    //     Ok(out)
    // }
    // COSMolKit✔️✔️: Original validation order and linear numeric conversion.

    validate_matrix(atom_count, &rows, "2D", &[2, 3], "2 or 3")?;
    let mut output = Vec::with_capacity(rows.len());
    // Pinned d892 extract_2d_coordinates, original 1e-12 inclusive threshold.
    // Complexity: linear conversion with one output allocation, matching source.
    for (row, values) in rows.into_iter().enumerate() {
        if values.len() == 3 {
            match params.z_policy {
                CoordinateZPolicy::Ignore => {}
                CoordinateZPolicy::RequireZero if values[2].abs() <= 1.0e-12 => {}
                CoordinateZPolicy::RequireZero => {
                    return Err(CoordinateInputError::NonZeroZ { row, z: values[2] });
                }
                CoordinateZPolicy::Error => {
                    return Err(CoordinateInputError::ThreeColumnsForbidden);
                }
            }
        }
        output.push([values[0], values[1]]);
    }
    Ok(output)
}
/// Interpret finite XYZ rows without changing is_3d or projecting dimensions.
pub fn coordinates_3d_from_input(
    atom_count: usize,
    rows: Vec<Vec<f64>>,
) -> Result<Vec<[f64; 3]>, CoordinateInputError> {
    // Pinned COSMolKit d892 original numeric ingress, verbatim source:
    // fn extract_3d_coordinates(
    //     coords: &Bound<'_, PyAny>,
    //     expected_rows: usize,
    // ) -> PyResult<Vec<[f64; 3]>> {
    //     extract_coordinate_matrix(coords, expected_rows, &[3], "3D coordinates").map(|rows| {
    //         rows.into_iter()
    //             .map(|row| [row[0], row[1], row[2]])
    //             .collect()
    //     })
    // }
    // COSMolKit✔️✔️: Original validation order and linear numeric conversion.

    validate_matrix(atom_count, &rows, "3D", &[3], "3")?;
    Ok(rows
        .into_iter()
        .map(|row| [row[0], row[1], row[2]])
        .collect())
}

fn source_dimension(block: &CoordinateBlock) -> Option<CoordinateDimension> {
    // Pinned COSMolKit d892 original Rust, verbatim source anchor:
    // fn source_coordinate_dim_for_block(
    //     coord_block: &CoordinateBlock,
    // ) -> Option<crate::CoordinateDimension> {
    //     if coord_block
    //         .conformers_3d
    //         .iter()
    //         .any(crate::Conformer3D::is_3d)
    //     {
    //         Some(crate::CoordinateDimension::ThreeD)
    //     } else if !coord_block.conformers_2d.is_empty() || !coord_block.conformers_3d.is_empty() {
    //         Some(crate::CoordinateDimension::TwoD)
    //     } else {
    //         None
    //     }
    // }
    // COSMolKit✔️✔️: Native source transition and scan/allocation shape reproduced.
    // Domain operates in the runtime-provided detached block; no whole-state clone.

    // Pinned d892 source_coordinate_dim_for_block, behavior and scan cost match.
    if block.conformers_3d.iter().any(Conformer3D::is_3d) {
        Some(CoordinateDimension::ThreeD)
    } else if !block.conformers_2d.is_empty() || !block.conformers_3d.is_empty() {
        Some(CoordinateDimension::TwoD)
    } else {
        None
    }
}

pub fn install_2d_coordinates(
    block: &mut CoordinateBlock,
    atom_count: usize,
    rows: Vec<Vec<f64>>,
    params: &Coordinate2DInputParams,
) -> Result<(), CoordinateInputError> {
    // Pinned COSMolKit d892 original Rust, verbatim source anchor:
    // pub(super) fn with_2d_coordinate_block_impl(coords: Vec<[f64; 2]>) -> Result<(), OperationError> {
    //     let atom_count = {
    //         let read = parts.begin_topology_read()?;
    //         read.num_atoms()
    //     };
    //     if coords.len() != atom_count {
    //         return Err(OperationError::InvalidInput {
    //             operation: &WITH_2D_COORDINATE_BLOCK_SPEC,
    //             message: "2D coordinate row count mismatch",
    //         });
    //     }
    //
    //     parts.with_coordinates_mut(|_parts, coord_block| {
    //         coord_block.conformers_2d.clear();
    //         coord_block
    //             .conformers_2d
    //             .push(crate::Conformer2D::new(0, coords));
    //         coord_block.source_coordinate_dim = Some(crate::CoordinateDimension::TwoD);
    //         Ok(())
    //     })?;
    //     parts.clear_cache(DerivedState::DRAWING);
    //     Ok(())
    // }
    // COSMolKit✔️✔️: Native source transition and scan/allocation shape reproduced.
    // Domain operates in the runtime-provided detached block; no whole-state clone.

    let coords = coordinates_2d_from_input(atom_count, rows, params)?;
    block.clear_2d_conformers();
    block.record_source_conformer_append(CoordinateDimension::TwoD)?;
    block.conformers_2d.push(Conformer2D::new(0, coords));
    block.source_coordinate_dim = Some(CoordinateDimension::TwoD);
    Ok(())
}
pub fn replace_3d_coordinates(
    block: &mut CoordinateBlock,
    atom_count: usize,
    rows: Vec<Vec<f64>>,
    params: &Replace3DCoordinatesParams,
) -> Result<(), CoordinateInputError> {
    // Pinned COSMolKit d892 original Rust, verbatim source anchor:
    // pub(super) fn with_3d_coordinates_impl(
    //     coords: Vec<[f64; 3]>,
    //     conformer_index: usize,
    // ) -> Result<(), OperationError> {
    //     let atom_count = {
    //         let read = parts.begin_topology_read()?;
    //         read.num_atoms()
    //     };
    //     if coords.len() != atom_count {
    //         return Err(OperationError::InvalidInput {
    //             operation: &WITH_3D_COORDINATES_SPEC,
    //             message: "3D conformer row count mismatch",
    //         });
    //     }
    //
    //     parts.with_coordinates_mut(|_parts, coord_block| {
    //         if conformer_index >= coord_block.conformers_3d.len() {
    //             return Err(OperationError::InvalidInput {
    //                 operation: &WITH_3D_COORDINATES_SPEC,
    //                 message: "3D conformer index out of range",
    //             });
    //         }
    //         let existing = &coord_block.conformers_3d[conformer_index];
    //         coord_block.conformers_3d[conformer_index] =
    //             crate::Conformer3D::new(existing.id(), coords, existing.is_3d());
    //         coord_block.source_coordinate_dim = source_coordinate_dim_for_block(coord_block);
    //         Ok(())
    //     })?;
    //     parts.clear_cache(DerivedState::DRAWING);
    //     Ok(())
    // }
    // COSMolKit❗❌: ROOT CK-3795 approved actual ThreeD ID selection instead of
    // historical vector position. Default0 selects ID0; a missing ID fails.
    // ID lookup scans O(conformers), replacing source O(1) positional lookup.
    // Existing is_3d, selected property reset, provenance and other sets follow source.
    // Domain operates in the runtime-provided detached block; no whole-state clone.

    let coords = coordinates_3d_from_input(atom_count, rows)?;
    let index = block
        .conformers_3d
        .iter()
        .position(|conformer| conformer.id() == params.conformer_id)
        .ok_or(CoordinateInputError::ConformerNotFound {
            conformer_id: params.conformer_id,
            count: block.conformers_3d.len(),
        })?;
    let existing = &block.conformers_3d[index];
    block.conformers_3d[index] = Conformer3D::new(existing.id(), coords, existing.is_3d());
    block.source_coordinate_dim = source_dimension(block);
    Ok(())
}
pub fn append_3d_conformer(
    block: &mut CoordinateBlock,
    atom_count: usize,
    rows: Vec<Vec<f64>>,
    params: &Coordinate3DInputParams,
) -> Result<usize, CoordinateInputError> {
    // Pinned COSMolKit d892 original Rust, verbatim source anchor:
    // pub(super) fn with_added_3d_conformer_impl(
    //     coords: Vec<[f64; 3]>,
    //     is_3d: bool,
    // ) -> Result<(), OperationError> {
    //     let atom_count = {
    //         let read = parts.begin_topology_read()?;
    //         read.num_atoms()
    //     };
    //     if coords.len() != atom_count {
    //         return Err(OperationError::InvalidInput {
    //             operation: &WITH_ADDED_3D_CONFORMER_SPEC,
    //             message: "3D conformer row count mismatch",
    //         });
    //     }
    //
    //     parts.with_coordinates_mut(|_parts, coord_block| {
    //         let next_id = coord_block
    //             .conformers_3d
    //             .iter()
    //             .map(crate::Conformer3D::id)
    //             .max()
    //             .map_or(0, |max_id| max_id + 1);
    //         coord_block
    //             .conformers_3d
    //             .push(crate::Conformer3D::new(next_id, coords, is_3d));
    //         coord_block.source_coordinate_dim = source_coordinate_dim_for_block(coord_block);
    //         Ok(())
    //     })?;
    //     parts.clear_cache(DerivedState::DRAWING);
    //     Ok(())
    // }
    // COSMolKit✔️✔️: Native source transition and scan/allocation shape reproduced.
    // Domain operates in the runtime-provided detached block; no whole-state clone.

    let coords = coordinates_3d_from_input(atom_count, rows)?;
    let id = match block.conformers_3d.iter().map(Conformer3D::id).max() {
        Some(max_id) => max_id
            .checked_add(1)
            .ok_or(CoordinateInputError::ConformerIdOverflow { max_id })?,
        None => 0,
    };
    let position = block.conformers_3d.len();
    block.record_source_conformer_append(CoordinateDimension::ThreeD)?;
    block
        .conformers_3d
        .push(Conformer3D::new(id, coords, params.is_3d));
    block.source_coordinate_dim = source_dimension(block);
    // The original Python report is vector position, even for noncontiguous IDs.
    Ok(position)
}
pub fn install_only_3d_conformer(
    block: &mut CoordinateBlock,
    atom_count: usize,
    rows: Vec<Vec<f64>>,
    params: &Coordinate3DInputParams,
) -> Result<usize, CoordinateInputError> {
    // Pinned COSMolKit d892 original Rust, verbatim source anchor:
    // pub(super) fn with_only_3d_conformer_impl(
    //     coords: Vec<[f64; 3]>,
    //     is_3d: bool,
    // ) -> Result<(), OperationError> {
    //     let atom_count = {
    //         let read = parts.begin_topology_read()?;
    //         read.num_atoms()
    //     };
    //     if coords.len() != atom_count {
    //         return Err(OperationError::InvalidInput {
    //             operation: &WITH_ONLY_3D_CONFORMER_SPEC,
    //             message: "3D conformer row count mismatch",
    //         });
    //     }
    //
    //     parts.with_coordinates_mut(|_parts, coord_block| {
    //         coord_block.conformers_3d.clear();
    //         coord_block
    //             .conformers_3d
    //             .push(crate::Conformer3D::new(0, coords, is_3d));
    //         coord_block.source_coordinate_dim = source_coordinate_dim_for_block(coord_block);
    //         Ok(())
    //     })?;
    //     parts.clear_cache(DerivedState::DRAWING);
    //     Ok(())
    // }
    // COSMolKit✔️✔️: Native source transition and scan/allocation shape reproduced.
    // Domain operates in the runtime-provided detached block; no whole-state clone.

    let coords = coordinates_3d_from_input(atom_count, rows)?;
    block.clear_3d_conformers();
    block.record_source_conformer_append(CoordinateDimension::ThreeD)?;
    block
        .conformers_3d
        .push(Conformer3D::new(0, coords, params.is_3d));
    block.source_coordinate_dim = source_dimension(block);
    Ok(0)
}
pub fn clear_3d_conformers(block: &mut CoordinateBlock) {
    // Pinned COSMolKit d892 original Rust, verbatim source anchor:
    // pub(super) fn with_cleared_3d_conformers_impl() -> Result<(), OperationError> {
    //     parts.with_coordinates_mut(|_parts, coord_block| {
    //         coord_block.conformers_3d.clear();
    //         coord_block.source_coordinate_dim = source_coordinate_dim_for_block(coord_block);
    //         Ok(())
    //     })?;
    //     parts.clear_cache(DerivedState::DRAWING);
    //     Ok(())
    // }
    // COSMolKit✔️✔️: Native source transition and scan/allocation shape reproduced.
    // Domain operates in the runtime-provided detached block; no whole-state clone.

    block.clear_3d_conformers();
    block.source_coordinate_dim = source_dimension(block);
}

#[derive(Clone, Debug, PartialEq, Eq, thiserror::Error)]
pub enum Coordinate3DReadError {
    #[error("no 3D conformer present with ID {conformer_id} ({count} stored)")]
    ConformerNotFound { conformer_id: usize, count: usize },
}

/// Borrow XYZ rows by dimension-scoped identity; never generate or select XY.
pub fn coordinates_3d_for_id(
    block: &CoordinateBlock,
    conformer_id: usize,
) -> Result<&[[f64; 3]], Coordinate3DReadError> {
    //     fn coordinates_3d<'py>(
    //         &self,
    //         py: Python<'py>,
    //         conformer_index: usize,
    //     ) -> PyResult<Bound<'py, PyAny>> {
    //         let Some(coords) = self.inner.conformers_3d().get(conformer_index) else {
    //             return Err(PyValueError::new_err(format!(
    //                 "no 3D conformer present at index {conformer_index}"
    //             )));
    //         };
    //         let rows: Vec<Vec<f64>> = coords
    //             .coordinates()
    //             .iter()
    //             .map(|p| vec![p[0], p[1], p[2]])
    //             .collect();
    //         PyArray2::from_vec2(py, &rows)
    //             .map(|array| array.into_any())
    //             .map_err(|err| PyValueError::new_err(format!("Molecule.coordinates_3d failed: {err}")))
    //     }
    // COSMolKit❗❌: Approved coordinate_selection_contract requires IDs instead
    // of historical positions. Explicit ID0 is the documented Python default.
    // This is a selection difference, not source-exact positional behavior.
    // Linear ID lookup replaces constant-time positional access; borrowing
    // avoids the original owner's copied row matrix. Binding alone copies rows.
    block
        .conformers_3d
        .iter()
        .find(|value| value.id() == conformer_id)
        .map(Conformer3D::coordinates)
        .ok_or(Coordinate3DReadError::ConformerNotFound {
            conformer_id,
            count: block.conformers_3d.len(),
        })
}
