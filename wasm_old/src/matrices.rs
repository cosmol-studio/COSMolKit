//! Typed read-only projection of canonical distance matrices.
use crate::Molecule;
use cosmolkit as ck;
impl Molecule {
    pub fn distance_matrix(&self) -> Result<ck::DenseMatrix, ck::MatrixError> {
        // COSMolKit❗✔️: self.inner.distance_matrix()
        self.inner.borrow().distance_matrix()
    }
    pub fn distance_matrix_with_params(
        &self,
        params: &ck::DistanceMatrixParams,
    ) -> Result<ck::DenseMatrix, ck::MatrixError> {
        // COSMolKit❗✔️: self.inner.distance_matrix_with_params(params)
        self.inner.borrow().distance_matrix_with_params(params)
    }
    pub fn distance_matrix_3d(&self) -> Result<ck::DenseMatrix, ck::MatrixError> {
        // COSMolKit❗✔️: self.inner.distance_matrix_3d()
        self.inner.borrow().distance_matrix_3d()
    }
    pub fn distance_matrix_3d_with_params(
        &self,
        params: &ck::DistanceMatrix3dParams,
    ) -> Result<ck::DenseMatrix, ck::MatrixError> {
        // COSMolKit❗✔️: self.inner.distance_matrix_3d_with_params(params)
        self.inner.borrow().distance_matrix_3d_with_params(params)
    }
}
#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn all_distance_policies_and_3d_selectors_match_owner_without_mutation() {
        for text in ["", "CCO", "C=O", "c1ccccc1", "C.O", "*"] {
            let m = Molecule::from_smiles(text).unwrap();
            let owner = ck::Molecule::from_smiles(text).unwrap();
            let before = m.inner.borrow().clone();
            assert_eq!(m.distance_matrix(), owner.distance_matrix());
            for use_bond_order in [false, true] {
                for use_atom_weights in [false, true] {
                    let p = ck::DistanceMatrixParams {
                        use_bond_order,
                        use_atom_weights,
                    };
                    assert_eq!(
                        m.distance_matrix_with_params(&p),
                        owner.distance_matrix_with_params(&p)
                    );
                }
            }
            assert_eq!(*m.inner.borrow(), before);
            assert!(matches!(
                m.distance_matrix_3d(),
                Err(ck::MatrixError::No3dConformer)
            ));
            assert!(matches!(
                m.distance_matrix_3d_with_params(&ck::DistanceMatrix3dParams {
                    conformer_id: Some(99),
                    use_atom_weights: false
                }),
                Err(ck::MatrixError::ConformerNotFound { conformer_id: 99 })
            ));
        }
        let block = "3\nmatrix\nC 0 0 0\nN 3 0 0\nO 3 4 0\n";
        let m = Molecule::from_xyz_block(block).unwrap();
        let owner = ck::Molecule::from_xyz_block(block).unwrap();
        let before = m.inner.borrow().clone();
        assert_eq!(m.distance_matrix_3d(), owner.distance_matrix_3d());
        for conformer_id in [None, Some(0), Some(99)] {
            for use_atom_weights in [false, true] {
                let p = ck::DistanceMatrix3dParams {
                    conformer_id,
                    use_atom_weights,
                };
                assert_eq!(
                    m.distance_matrix_3d_with_params(&p),
                    owner.distance_matrix_3d_with_params(&p)
                );
            }
        }
        let matrix = m.distance_matrix_3d().unwrap();
        assert_eq!(matrix.dimension(), 3);
        assert_eq!(matrix.values(), [0., 3., 5., 3., 0., 4., 5., 4., 0.]);
        assert_eq!(matrix.get(3, 0), None);
        assert_eq!(matrix.get(usize::MAX, usize::MAX), None);
        assert_eq!(*m.inner.borrow(), before);
    }
}
