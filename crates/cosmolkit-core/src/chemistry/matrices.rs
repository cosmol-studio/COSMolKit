//! Shared source-backed molecular distance matrices.

use std::sync::Arc;

use crate::Molecule;

pub(crate) const LOCAL_INF_DISTANCE: f64 = 1.0e8;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum DistanceMatrixError {
    BondOrderWeightingOutsideSelectedBoundary,
    AtomWeightingOutsideSelectedBoundary,
    MissingConformer {
        conf_id: i32,
    },
    CoordinateCountMismatch {
        conf_id: i32,
        expected: usize,
        actual: usize,
    },
}

/// Reproduce the selected unweighted `MolOps::getDistanceMat()` boundary.
pub(crate) fn distance_matrix(
    molecule: &Molecule,
    use_bond_order: bool,
    use_atom_weights: bool,
) -> Result<Arc<[f64]>, DistanceMatrixError> {
    // RDKit✔️✔️: double *getDistanceMat(const ROMol &mol, bool useBO, bool useAtomWts,
    // RDKit✔️✔️:                        bool force, const char *propNamePrefix) {
    // RDKit✔️✔️:   std::string propName;
    // RDKit✔️✔️:   boost::shared_array<double> sptr;
    // RDKit✔️✔️:   if (propNamePrefix) {
    // RDKit✔️✔️:     propName = propNamePrefix;
    // RDKit✔️✔️:   } else {
    // RDKit✔️✔️:     propName = "";
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   propName += "DistanceMatrix";
    // RDKit✔️✔️:   // make sure we don't use the nonBO cache for the BO matrix and vice versa:
    // RDKit✔️✔️:   if (useBO) {
    // RDKit✔️✔️:     propName += "BO";
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   if (useAtomWts) {
    // RDKit✔️✔️:     propName += "AtomWts";
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   if (!force && mol.hasProp(propName)) {
    // RDKit✔️✔️:     mol.getProp(propName, sptr);
    // RDKit✔️✔️:     return sptr.get();
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   int nAts = mol.getNumAtoms();
    // RDKit✔️✔️:   auto *dMat = new double[nAts * nAts];
    // RDKit✔️✔️:   int i, j;
    // RDKit✔️✔️:   // initialize off diagonals to LOCAL_INF and diagonals to 0
    // RDKit✔️✔️:   for (i = 0; i < nAts * nAts; i++) {
    // RDKit✔️✔️:     dMat[i] = LOCAL_INF;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   for (i = 0; i < nAts; i++) {
    // RDKit✔️✔️:     dMat[i * nAts + i] = 0.0;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   ROMol::EDGE_ITER firstB, lastB;
    // RDKit✔️✔️:   boost::tie(firstB, lastB) = mol.getEdges();
    // RDKit✔️✔️:   while (firstB != lastB) {
    // RDKit✔️✔️:     const Bond *bond = mol[*firstB];
    // RDKit✔️✔️:     i = bond->getBeginAtomIdx();
    // RDKit✔️✔️:     j = bond->getEndAtomIdx();
    // RDKit✔️✔️:     double contrib;
    // RDKit❌❌:     if (useBO) {
    // RDKit❌❌:       if (!bond->getIsAromatic()) {
    // RDKit❌❌:         contrib = 1. / bond->getBondTypeAsDouble();
    // RDKit❌❌:       } else {
    // RDKit❌❌:         contrib = 2. / 3.;
    // RDKit❌❌:       }
    // RDKit✔️✔️:     } else {
    // RDKit✔️✔️:       contrib = 1.0;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     dMat[i * nAts + j] = contrib;
    // RDKit✔️✔️:     dMat[j * nAts + i] = contrib;
    // RDKit✔️✔️:     ++firstB;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   auto *pathMat = new int[nAts * nAts];
    // RDKit✔️✔️:   memset(static_cast<void *>(pathMat), 0, nAts * nAts * sizeof(int));
    // RDKit✔️✔️:   FloydWarshall(nAts, dMat, pathMat);
    // RDKit✔️✔️:
    // RDKit❌❌:   if (useAtomWts) {
    // RDKit❌❌:     for (i = 0; i < nAts; i++) {
    // RDKit❌❌:       int anum = mol.getAtomWithIdx(i)->getAtomicNum();
    // RDKit❌❌:       dMat[i * nAts + i] = 6.0 / anum;
    // RDKit❌❌:     }
    // RDKit❌❌:   }
    // RDKit✔️✔️:   sptr.reset(dMat);
    // RDKit✔️✔️:   mol.setProp(propName, sptr, true);
    // RDKit✔️🔝:   boost::shared_array<int> iSptr(pathMat);
    // RDKit✔️🔝:   mol.setProp(propName + "_Paths", iSptr, true);
    // RDKit✔️✔️:
    // RDKit✔️✔️:   return dMat;
    // RDKit✔️✔️: };
    //
    // The source path matrix is an internal computed property that no modeled
    // COSMolKit API reads. The Floyd-Warshall port still computes the same path
    // entries while producing the distance matrix, then releases that second
    // O(n^2) buffer instead of retaining unobservable state. The distance matrix
    // itself is retained in the molecule's logically read-only computed cache.
    if use_bond_order {
        return Err(DistanceMatrixError::BondOrderWeightingOutsideSelectedBoundary);
    }
    if use_atom_weights {
        return Err(DistanceMatrixError::AtomWeightingOutsideSelectedBoundary);
    }

    Ok(molecule.topological_distance_matrix_cache_or_init(|| {
        let dimension = molecule.num_atoms();
        let mut last_distances = vec![LOCAL_INF_DISTANCE; dimension * dimension];
        let mut last_paths = vec![0_i32; dimension * dimension];

        for atom_index in 0..dimension {
            last_distances[atom_index * dimension + atom_index] = 0.0;
        }
        for bond in molecule.bonds() {
            let begin = bond.begin().index();
            let end = bond.end().index();
            last_distances[begin * dimension + end] = 1.0;
            last_distances[end * dimension + begin] = 1.0;
        }

        for row in 0..dimension {
            let row_offset = row * dimension;
            for column in 0..dimension {
                last_paths[row_offset + column] =
                    if row == column || last_distances[row_offset + column] == LOCAL_INF_DISTANCE {
                        -1
                    } else {
                        row as i32
                    };
            }
        }

        let mut current_distances = vec![0.0; dimension * dimension];
        let mut current_paths = vec![0_i32; dimension * dimension];
        // RDKit✔️✔️: template <class T>
        // RDKit✔️✔️: void FloydWarshall(int dim, T *adjMat, int *pathMat) {
        // RDKit✔️✔️:   int k, i, j;
        // RDKit✔️✔️:   T *currD, *lastD, *tTemp;
        // RDKit✔️✔️:   int *currP, *lastP, *iTemp;
        // RDKit✔️✔️:
        // RDKit✔️✔️:   currD = new T[dim * dim];
        // RDKit✔️✔️:   currP = new int[dim * dim];
        // RDKit✔️✔️:   lastD = new T[dim * dim];
        // RDKit✔️✔️:   lastP = new int[dim * dim];
        // RDKit✔️✔️:
        // RDKit✔️✔️:   memcpy(static_cast<void *>(lastD), static_cast<void *>(adjMat),
        // RDKit✔️✔️:          dim * dim * sizeof(T));
        // RDKit✔️✔️:
        // RDKit✔️✔️:   // initialize the paths
        // RDKit✔️✔️:   for (i = 0; i < dim; i++) {
        // RDKit✔️✔️:     int itab = i * dim;
        // RDKit✔️✔️:     for (j = 0; j < dim; j++) {
        // RDKit✔️✔️:       if (i == j || adjMat[itab + j] == LOCAL_INF) {
        // RDKit✔️✔️:         pathMat[itab + j] = -1;
        // RDKit✔️✔️:       } else {
        // RDKit✔️✔️:         pathMat[itab + j] = i;
        // RDKit✔️✔️:       }
        // RDKit✔️✔️:     }
        // RDKit✔️✔️:   }
        // RDKit✔️✔️:   memcpy(static_cast<void *>(lastP), static_cast<void *>(pathMat),
        // RDKit✔️✔️:          dim * dim * sizeof(int));
        // RDKit✔️✔️:
        // RDKit✔️✔️:   for (k = 0; k < dim; k++) {
        // RDKit✔️✔️:     int ktab = k * dim;
        // RDKit✔️✔️:     for (i = 0; i < dim; i++) {
        // RDKit✔️✔️:       int itab = i * dim;
        // RDKit✔️✔️:       for (j = 0; j < dim; j++) {
        // RDKit✔️✔️:         T v1 = lastD[itab + j];
        // RDKit✔️✔️:         T v2 = lastD[itab + k] + lastD[ktab + j];
        // RDKit✔️✔️:         if (v1 <= v2) {
        // RDKit✔️✔️:           currD[itab + j] = v1;
        // RDKit✔️✔️:           currP[itab + j] = lastP[itab + j];
        // RDKit✔️✔️:         } else {
        // RDKit✔️✔️:           currD[itab + j] = v2;
        // RDKit✔️✔️:           currP[itab + j] = lastP[ktab + j];
        // RDKit✔️✔️:         }
        // RDKit✔️✔️:       }
        // RDKit✔️✔️:     }
        // RDKit✔️✔️:     tTemp = currD;
        // RDKit✔️✔️:     currD = lastD;
        // RDKit✔️✔️:     lastD = tTemp;
        // RDKit✔️✔️:
        // RDKit✔️✔️:     iTemp = currP;
        // RDKit✔️✔️:     currP = lastP;
        // RDKit✔️✔️:     lastP = iTemp;
        // RDKit✔️✔️:   }
        // RDKit✔️✔️:   memcpy(static_cast<void *>(adjMat), static_cast<void *>(lastD),
        // RDKit✔️✔️:          dim * dim * sizeof(T));
        // RDKit✔️✔️:   memcpy(static_cast<void *>(pathMat), static_cast<void *>(lastP),
        // RDKit✔️✔️:          dim * dim * sizeof(int));
        // RDKit✔️✔️:
        // RDKit✔️✔️:   delete[] currD;
        // RDKit✔️✔️:   delete[] currP;
        // RDKit✔️✔️:   delete[] lastD;
        // RDKit✔️✔️:   delete[] lastP;
        // RDKit✔️✔️: }
        for intermediate in 0..dimension {
            let intermediate_offset = intermediate * dimension;
            for row in 0..dimension {
                let row_offset = row * dimension;
                for column in 0..dimension {
                    let direct = last_distances[row_offset + column];
                    let through =
                        last_distances[row_offset + intermediate] + last_distances[intermediate_offset + column];
                    if direct <= through {
                        current_distances[row_offset + column] = direct;
                        current_paths[row_offset + column] = last_paths[row_offset + column];
                    } else {
                        current_distances[row_offset + column] = through;
                        current_paths[row_offset + column] = last_paths[intermediate_offset + column];
                    }
                }
            }
            std::mem::swap(&mut current_distances, &mut last_distances);
            std::mem::swap(&mut current_paths, &mut last_paths);
        }

        last_distances
    }))
}

#[must_use]
pub(crate) fn topological_distance_matrix(molecule: &Molecule) -> Arc<[f64]> {
    distance_matrix(molecule, false, false).expect("the unweighted distance-matrix boundary is infallible")
}

/// Reproduce `MolOps::get3DDistanceMat()` for the selected unweighted boundary.
pub(crate) fn distance_matrix_3d(
    molecule: &Molecule,
    conf_id: i32,
    use_atom_weights: bool,
) -> Result<Arc<[f64]>, DistanceMatrixError> {
    // RDKit✔️✔️: double *get3DDistanceMat(const ROMol &mol, int confId, bool useAtomWts,
    // RDKit✔️✔️:                          bool force, const char *propNamePrefix) {
    // RDKit✔️✔️:   const Conformer &conf = mol.getConformer(confId);
    // RDKit✔️✔️:   std::string propName;
    // RDKit✔️✔️:   boost::shared_array<double> sptr;
    // RDKit✔️✔️:   if (propNamePrefix) {
    // RDKit✔️✔️:     propName = propNamePrefix;
    // RDKit✔️✔️:   } else {
    // RDKit✔️✔️:     propName = "_";
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   if (propName != "") {
    // RDKit✔️✔️:     propName += "3DDistanceMatrix_Conf" + std::to_string(conf.getId());
    // RDKit❌❌:     if (useAtomWts) {
    // RDKit❌❌:       propName += "_AtomWeights";
    // RDKit❌❌:     }
    // RDKit✔️✔️:     if (!force && mol.hasProp(propName)) {
    // RDKit✔️✔️:       mol.getProp(propName, sptr);
    // RDKit✔️✔️:       return sptr.get();
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   unsigned int nAts = mol.getNumAtoms();
    // RDKit✔️✔️:   auto *dMat = new double[nAts * nAts];
    // RDKit✔️✔️:
    // RDKit✔️✔️:   for (unsigned int i = 0; i < nAts; ++i) {
    // RDKit❌❌:     if (useAtomWts) {
    // RDKit❌❌:       dMat[i * nAts + i] = 6.0 / mol.getAtomWithIdx(i)->getAtomicNum();
    // RDKit✔️✔️:     } else {
    // RDKit✔️✔️:       dMat[i * nAts + i] = 0.0;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     for (unsigned int j = i + 1; j < nAts; ++j) {
    // RDKit✔️✔️:       double dist = (conf.getAtomPos(i) - conf.getAtomPos(j)).length();
    // RDKit✔️✔️:       dMat[i * nAts + j] = dist;
    // RDKit✔️✔️:       dMat[j * nAts + i] = dist;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   if (propName != "") {
    // RDKit✔️✔️:     sptr.reset(dMat);
    // RDKit✔️✔️:     mol.setProp(propName, sptr, true);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return dMat;
    // RDKit✔️✔️: }
    //
    // RDKit✔️✔️: const Conformer &ROMol::getConformer(int id) const {
    // RDKit✔️✔️:   // make sure we have more than one conformation
    // RDKit✔️✔️:   if (d_confs.size() == 0) {
    // RDKit✔️✔️:     throw ConformerException("No conformations available on the molecule");
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   if (id < 0) {
    // RDKit✔️✔️:     return *(d_confs.front());
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   auto cid = (unsigned int)id;
    // RDKit✔️✔️:   for (auto conf : d_confs) {
    // RDKit✔️✔️:     if (conf->getId() == cid) {
    // RDKit✔️✔️:       return *conf;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   // we did not find a conformation with the specified ID
    // RDKit✔️✔️:   std::string mesg = "Can't find conformation with ID: ";
    // RDKit✔️✔️:   mesg += id;
    // RDKit✔️✔️:   throw ConformerException(mesg);
    // RDKit✔️✔️: }
    if use_atom_weights {
        return Err(DistanceMatrixError::AtomWeightingOutsideSelectedBoundary);
    }

    let conformer = if conf_id < 0 {
        molecule.conformers_3d().first()
    } else {
        let requested = conf_id as usize;
        molecule
            .conformers_3d()
            .iter()
            .find(|conformer| conformer.id() == requested)
    }
    .ok_or(DistanceMatrixError::MissingConformer { conf_id })?;

    if conformer.coordinates().len() != molecule.num_atoms() {
        return Err(DistanceMatrixError::CoordinateCountMismatch {
            conf_id,
            expected: molecule.num_atoms(),
            actual: conformer.coordinates().len(),
        });
    }
    Ok(molecule.distance_matrix_3d_cache_or_init(conformer.id(), || {
        distance_matrix_3d_from_coordinates(conformer.coordinates(), molecule.num_atoms(), conf_id)
            .expect("coordinate dimensions were validated before caching")
    }))
}

fn distance_matrix_3d_from_coordinates(
    coordinates: &[[f64; 3]],
    atom_count: usize,
    conf_id: i32,
) -> Result<Vec<f64>, DistanceMatrixError> {
    if coordinates.len() != atom_count {
        return Err(DistanceMatrixError::CoordinateCountMismatch {
            conf_id,
            expected: atom_count,
            actual: coordinates.len(),
        });
    }

    let mut distances = vec![0.0; atom_count * atom_count];
    for left in 0..atom_count {
        for right in left + 1..atom_count {
            let delta_x = coordinates[left][0] - coordinates[right][0];
            let delta_y = coordinates[left][1] - coordinates[right][1];
            let delta_z = coordinates[left][2] - coordinates[right][2];
            let distance = (delta_x * delta_x + delta_y * delta_y + delta_z * delta_z).sqrt();
            distances[left * atom_count + right] = distance;
            distances[right * atom_count + left] = distance;
        }
    }
    Ok(distances)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{AtomPairFingerprintParams, AtomSpec, Conformer3D, Element};

    fn two_atom_molecule_with_conformers(conformers: impl IntoIterator<Item = Conformer3D>) -> Molecule {
        let mut builder = Molecule::builder();
        builder.add_atom(AtomSpec::new(Element::C));
        builder.add_atom(AtomSpec::new(Element::O));
        for conformer in conformers {
            builder.add_conformer(conformer).unwrap();
        }
        builder.build().unwrap()
    }

    #[test]
    fn rdkit_2d_distance_empty_molecule_is_an_empty_row_major_matrix() {
        assert!(topological_distance_matrix(&Molecule::new()).is_empty());
    }

    #[test]
    fn rdkit_2d_distance_preserves_source_atom_and_row_order() {
        let molecule = Molecule::from_smiles_with_sanitize("CCC", false).unwrap();
        assert_eq!(
            topological_distance_matrix(&molecule).as_ref(),
            &[0.0, 1.0, 2.0, 1.0, 0.0, 1.0, 2.0, 1.0, 0.0]
        );
    }

    #[test]
    fn rdkit_2d_distance_keeps_disconnected_pairs_at_local_inf() {
        let molecule = Molecule::from_smiles_with_sanitize("CC.C", false).unwrap();
        assert_eq!(
            topological_distance_matrix(&molecule).as_ref(),
            &[
                0.0,
                1.0,
                LOCAL_INF_DISTANCE,
                1.0,
                0.0,
                LOCAL_INF_DISTANCE,
                LOCAL_INF_DISTANCE,
                LOCAL_INF_DISTANCE,
                0.0,
            ]
        );
    }

    #[test]
    fn rdkit_2d_distance_rejects_unselected_weighted_options() {
        let molecule = Molecule::from_smiles_with_sanitize("CO", false).unwrap();
        assert_eq!(
            distance_matrix(&molecule, true, false),
            Err(DistanceMatrixError::BondOrderWeightingOutsideSelectedBoundary)
        );
        assert_eq!(
            distance_matrix(&molecule, false, true),
            Err(DistanceMatrixError::AtomWeightingOutsideSelectedBoundary)
        );
        assert_eq!(
            distance_matrix(&molecule, true, true),
            Err(DistanceMatrixError::BondOrderWeightingOutsideSelectedBoundary)
        );
    }

    #[test]
    fn distgeom_uses_the_shared_distance_matrix_without_a_local_copy() {
        let source = include_str!("distgeom.rs");
        assert!(source.contains("topological_distance_matrix as compute_topological_distances"));
        assert_eq!(source.matches("fn compute_topological_distances(").count(), 0);
        assert_eq!(source.matches("const LOCAL_INF_DIST").count(), 0);
    }

    #[test]
    fn rdkit_3d_distance_default_and_any_negative_id_select_first_conformer() {
        let molecule = two_atom_molecule_with_conformers([
            Conformer3D::new(7, vec![[0.0, 0.0, 0.0], [3.0, 4.0, 0.0]], true),
            Conformer3D::new(2, vec![[0.0, 0.0, 0.0], [0.0, 0.0, 2.0]], true),
        ]);
        let expected = vec![0.0, 5.0, 5.0, 0.0];
        assert_eq!(distance_matrix_3d(&molecule, -1, false).unwrap().as_ref(), expected);
        assert_eq!(distance_matrix_3d(&molecule, -9, false).unwrap().as_ref(), expected);
    }

    #[test]
    fn rdkit_3d_distance_explicit_id_selects_matching_conformer_not_position() {
        let molecule = two_atom_molecule_with_conformers([
            Conformer3D::new(7, vec![[0.0, 0.0, 0.0], [3.0, 4.0, 0.0]], true),
            Conformer3D::new(2, vec![[0.0, 0.0, 0.0], [0.0, 0.0, 2.0]], true),
        ]);
        assert_eq!(
            distance_matrix_3d(&molecule, 2, false).unwrap().as_ref(),
            &[0.0, 2.0, 2.0, 0.0]
        );
    }

    #[test]
    fn rdkit_3d_distance_empty_molecule_with_empty_conformer_is_empty() {
        let mut builder = Molecule::builder();
        builder.add_conformer(Conformer3D::new(4, Vec::new(), true)).unwrap();
        let molecule = builder.build().unwrap();
        assert!(distance_matrix_3d(&molecule, -1, false).unwrap().is_empty());
    }

    #[test]
    fn rdkit_3d_distance_preserves_fractional_euclidean_arithmetic_and_order() {
        let molecule =
            two_atom_molecule_with_conformers([Conformer3D::new(0, vec![[0.25, -0.5, 1.5], [1.75, 1.5, -0.5]], true)]);
        let expected = (1.5_f64 * 1.5 + 2.0 * 2.0 + 2.0 * 2.0).sqrt();
        assert_eq!(
            distance_matrix_3d(&molecule, 0, false).unwrap().as_ref(),
            &[0.0, expected, expected, 0.0]
        );
    }

    #[test]
    fn rdkit_3d_distance_reports_missing_conformers_and_unselected_weighting() {
        let mut builder = Molecule::builder();
        builder.add_atom(AtomSpec::new(Element::C));
        let no_conformer = builder.build().unwrap();
        assert_eq!(
            distance_matrix_3d(&no_conformer, -1, false),
            Err(DistanceMatrixError::MissingConformer { conf_id: -1 })
        );

        let molecule =
            two_atom_molecule_with_conformers([Conformer3D::new(7, vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]], true)]);
        assert_eq!(
            distance_matrix_3d(&molecule, 3, false),
            Err(DistanceMatrixError::MissingConformer { conf_id: 3 })
        );
        assert_eq!(
            distance_matrix_3d(&molecule, -1, true),
            Err(DistanceMatrixError::AtomWeightingOutsideSelectedBoundary)
        );
    }

    #[test]
    fn rdkit_3d_distance_rejects_malformed_coordinate_count_before_indexing() {
        assert_eq!(
            distance_matrix_3d_from_coordinates(&[[0.0, 0.0, 0.0]], 2, 7),
            Err(DistanceMatrixError::CoordinateCountMismatch {
                conf_id: 7,
                expected: 2,
                actual: 1,
            })
        );
    }

    #[test]
    fn rdkit_3d_distance_propagates_non_finite_coordinate_arithmetic_like_source() {
        let distances = distance_matrix_3d_from_coordinates(&[[0.0, 0.0, 0.0], [f64::NAN, 0.0, 0.0]], 2, 0).unwrap();
        assert!(distances[1].is_nan());
        assert!(distances[2].is_nan());
        assert_eq!(distances[0], 0.0);
        assert_eq!(distances[3], 0.0);
    }

    #[test]
    fn rdkit_2d_distance_reuses_one_logically_read_only_cache_entry() {
        let molecule = Molecule::from_smiles("C1CCCCCCCCCCCCCCC1").unwrap();
        let cold = topological_distance_matrix(&molecule);
        let warm = topological_distance_matrix(&molecule);

        assert!(Arc::ptr_eq(&cold, &warm));
    }

    #[test]
    fn rdkit_2d_distance_cache_is_clone_local_and_invalidated_by_topology_operations() {
        let molecule = Molecule::from_smiles("C").unwrap();
        let source_cached = topological_distance_matrix(&molecule);
        let clone = molecule.clone();
        let inherited = topological_distance_matrix(&clone);
        assert!(Arc::ptr_eq(&source_cached, &inherited));

        let hydrogenated = clone.with_hydrogens().unwrap();
        let changed = topological_distance_matrix(&hydrogenated);
        assert_eq!(changed.len(), hydrogenated.num_atoms() * hydrogenated.num_atoms());
        assert!(!Arc::ptr_eq(&inherited, &changed));

        let source_warm = topological_distance_matrix(&molecule);
        assert!(Arc::ptr_eq(&source_cached, &source_warm));
    }

    #[test]
    fn rdkit_3d_distance_cache_uses_resolved_conformer_id() {
        let molecule = two_atom_molecule_with_conformers([
            Conformer3D::new(7, vec![[0.0, 0.0, 0.0], [3.0, 4.0, 0.0]], true),
            Conformer3D::new(2, vec![[0.0, 0.0, 0.0], [0.0, 0.0, 2.0]], true),
        ]);

        let default = distance_matrix_3d(&molecule, -1, false).unwrap();
        let same_resolved_id = distance_matrix_3d(&molecule, 7, false).unwrap();
        let second = distance_matrix_3d(&molecule, 2, false).unwrap();
        let second_warm = distance_matrix_3d(&molecule, 2, false).unwrap();

        assert!(Arc::ptr_eq(&default, &same_resolved_id));
        assert!(!Arc::ptr_eq(&default, &second));
        assert!(Arc::ptr_eq(&second, &second_warm));
    }

    #[test]
    fn rdkit_3d_distance_cache_is_invalidated_by_coordinate_operations() {
        let molecule =
            two_atom_molecule_with_conformers([Conformer3D::new(0, vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]], true)]);
        let before = distance_matrix_3d(&molecule, 0, false).unwrap();
        let moved = molecule.with_atom_position(1, [0.0, 0.0, 6.0], 0).unwrap();
        let after = distance_matrix_3d(&moved, 0, false).unwrap();

        assert!(!Arc::ptr_eq(&before, &after));
        assert_eq!(after.as_ref(), &[0.0, 6.0, 6.0, 0.0]);
        assert!(Arc::ptr_eq(&before, &distance_matrix_3d(&molecule, 0, false).unwrap()));
    }

    #[test]
    fn rdkit_distance_cache_parallel_initialization_publishes_one_matrix() {
        let molecule = Arc::new(Molecule::from_smiles("C1CCCCCCCCCCCCCCC1").unwrap());
        let matrices = std::thread::scope(|scope| {
            let handles = (0..16)
                .map(|_| {
                    let molecule = Arc::clone(&molecule);
                    scope.spawn(move || topological_distance_matrix(&molecule))
                })
                .collect::<Vec<_>>();
            handles
                .into_iter()
                .map(|handle| handle.join().unwrap())
                .collect::<Vec<_>>()
        });

        assert!(matrices.iter().skip(1).all(|matrix| Arc::ptr_eq(&matrices[0], matrix)));
    }

    #[test]
    fn rdkit_atom_pair_exact_output_is_unchanged_between_cold_and_warm_calls() {
        let molecule = Molecule::from_smiles("CCO").unwrap();
        let params = AtomPairFingerprintParams::default();
        let cold = molecule.atom_pair_fingerprint(&params).unwrap();
        let warm = molecule.atom_pair_fingerprint(&params).unwrap();

        assert_eq!(cold.on_bits(), &[1404, 1596, 1916]);
        assert_eq!(warm, cold);
    }
}
