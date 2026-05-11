use crate::Molecule;

#[derive(Debug, Clone, PartialEq, Eq, thiserror::Error)]
pub enum DgBoundsError {
    #[error("distance-geometry bounds matrix is not implemented")]
    NotImplemented,
    #[error(transparent)]
    UnsupportedFeature(#[from] crate::UnsupportedFeatureError),
}

pub fn dg_bounds_matrix(_molecule: &Molecule) -> Result<Vec<Vec<f64>>, DgBoundsError> {
    Err(crate::UnsupportedFeatureError::from_spec(&crate::COORDINATE_2D_FEATURE).into())
}
