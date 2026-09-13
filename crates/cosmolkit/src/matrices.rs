//! Thin read-only matrix queries over the detached core owner.

use crate::Molecule;

/// Public options for a full-molecule topological distance matrix.
///
/// Active atom and bond subsets remain detached implementation APIs because
/// they do not describe an ordinary single-molecule query projection.
#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub struct DistanceMatrixParams {
    pub use_bond_order: bool,
    pub use_atom_weights: bool,
}

impl Molecule {
    /// Returns the full topological distance matrix with source defaults.
    pub fn distance_matrix(
        &self,
    ) -> Result<cosmolkit_core::DenseMatrix, cosmolkit_core::MatrixError> {
        self.distance_matrix_with_params(&DistanceMatrixParams::default())
    }

    /// Returns the full topological distance matrix with explicit weighting.
    pub fn distance_matrix_with_params(
        &self,
        params: &DistanceMatrixParams,
    ) -> Result<cosmolkit_core::DenseMatrix, cosmolkit_core::MatrixError> {
        let detached = cosmolkit_core::TopologicalDistanceMatrixParams {
            use_bond_order: params.use_bond_order,
            use_atom_weights: params.use_atom_weights,
            active_atoms: None,
            active_bonds: None,
        };
        cosmolkit_core::topological_distance_matrix(self.topology(), &detached)
    }

    /// Returns the 3D distance matrix with source defaults.
    pub fn distance_matrix_3d(
        &self,
    ) -> Result<cosmolkit_core::DenseMatrix, cosmolkit_core::MatrixError> {
        self.distance_matrix_3d_with_params(&cosmolkit_core::DistanceMatrix3dParams::default())
    }

    /// Returns the 3D distance matrix for the selected conformer and weighting.
    pub fn distance_matrix_3d_with_params(
        &self,
        params: &cosmolkit_core::DistanceMatrix3dParams,
    ) -> Result<cosmolkit_core::DenseMatrix, cosmolkit_core::MatrixError> {
        cosmolkit_core::distance_matrix_3d(self.topology(), self.coordinates(), params)
    }
}
