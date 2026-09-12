//! Detached 2D layout and depiction boundaries.

use cosmolkit_model::{CoordinateBlock, TopologyBlock};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DepictError {
    Unsupported,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct DepictOptions {
    pub width: u32,
    pub height: u32,
}

pub fn layout_2d(
    topology: &TopologyBlock,
    options: &DepictOptions,
) -> Result<CoordinateBlock, DepictError> {
    let _ = (topology, options);
    Err(DepictError::Unsupported)
}

pub fn render_svg(
    topology: &TopologyBlock,
    coordinates: &CoordinateBlock,
    options: &DepictOptions,
) -> Result<String, DepictError> {
    let _ = (topology, coordinates, options);
    Err(DepictError::Unsupported)
}
