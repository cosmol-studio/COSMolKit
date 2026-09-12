//! Source-backed dense matrices over detached topology and coordinate values.
//!
//! This module has no live-molecule or cache authority. Runtime cache keys and
//! invalidation remain the sole responsibility of `cosmolkit`.

use std::collections::BTreeMap;

use cosmolkit_model::{
    AtomId, Bond, BondId, BondOrder, CoordinateBlock, CoordinateValidationError, TopologyBlock,
    TopologyValidationError,
};

use crate::{bond_type_as_double, bond_valence_contrib};

const LOCAL_INF: f64 = 100_000_000.0;

#[derive(Debug, Clone, PartialEq)]
pub struct DenseMatrix {
    dimension: usize,
    values: Vec<f64>,
}

impl DenseMatrix {
    #[must_use]
    pub const fn dimension(&self) -> usize {
        self.dimension
    }

    #[must_use]
    pub fn values(&self) -> &[f64] {
        &self.values
    }

    #[must_use]
    pub fn get(&self, row: usize, column: usize) -> Option<f64> {
        if row >= self.dimension || column >= self.dimension {
            return None;
        }
        self.values.get(row * self.dimension + column).copied()
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct AdjacencyMatrixParams {
    pub use_bond_order: bool,
    pub empty_value: i32,
    pub bonds_to_use: Option<Vec<BondId>>,
}

impl Default for AdjacencyMatrixParams {
    fn default() -> Self {
        Self {
            use_bond_order: false,
            empty_value: 0,
            bonds_to_use: None,
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub struct TopologicalDistanceMatrixParams {
    pub use_bond_order: bool,
    pub use_atom_weights: bool,
    pub active_atoms: Option<Vec<AtomId>>,
    pub active_bonds: Option<Vec<BondId>>,
}

#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub struct DistanceMatrix3dParams {
    pub conformer_id: Option<usize>,
    pub use_atom_weights: bool,
}

#[derive(Debug, Clone, PartialEq, Eq, thiserror::Error)]
pub enum MatrixError {
    #[error(transparent)]
    InvalidTopology(TopologyValidationError),
    #[error(transparent)]
    InvalidCoordinates(CoordinateValidationError),
    #[error("active bonds require an active atom selection")]
    ActiveBondsWithoutAtoms,
    #[error("active atom at position {position} is {atom}, outside {atom_count} atoms")]
    ActiveAtomOutOfRange {
        position: usize,
        atom: AtomId,
        atom_count: usize,
    },
    #[error("active atom {atom} is repeated at positions {first_position} and {second_position}")]
    DuplicateActiveAtom {
        atom: AtomId,
        first_position: usize,
        second_position: usize,
    },
    #[error("active bond at position {position} is {bond}, outside {bond_count} bonds")]
    ActiveBondOutOfRange {
        position: usize,
        bond: BondId,
        bond_count: usize,
    },
    #[error("active bond {bond} is repeated at positions {first_position} and {second_position}")]
    DuplicateActiveBond {
        bond: BondId,
        first_position: usize,
        second_position: usize,
    },
    #[error("active bond {bond} {endpoint} endpoint atom {atom} is absent from active atoms")]
    ActiveBondEndpointMissing {
        bond: BondId,
        endpoint: &'static str,
        atom: AtomId,
    },
    #[error("bond {bond} has source-unsupported matrix bond order {order:?}")]
    UnsupportedBondOrder { bond: BondId, order: BondOrder },
    #[error("no 3D conformer is available")]
    No3dConformer,
    #[error("3D conformer id {conformer_id} was not found")]
    ConformerNotFound { conformer_id: usize },
    #[error("matrix dimension {dimension} overflows row-major storage")]
    MatrixDimensionOverflow { dimension: usize },
}

fn matrix_len(dimension: usize) -> Result<usize, MatrixError> {
    dimension
        .checked_mul(dimension)
        .ok_or(MatrixError::MatrixDimensionOverflow { dimension })
}

fn source_bond_order(bond: &Bond) -> Result<f64, MatrixError> {
    bond_type_as_double(bond.order()).map_err(|_| MatrixError::UnsupportedBondOrder {
        bond: bond.id(),
        order: bond.order(),
    })
}

fn source_valence_contribution(bond: &Bond, atom: AtomId) -> Result<f64, MatrixError> {
    bond_valence_contrib(bond, atom).map_err(|_| MatrixError::UnsupportedBondOrder {
        bond: bond.id(),
        order: bond.order(),
    })
}

fn validated_bond_filter(
    topology: &TopologyBlock,
    selected: Option<&[BondId]>,
) -> Result<Option<Vec<bool>>, MatrixError> {
    let Some(selected) = selected else {
        return Ok(None);
    };
    let mut included = vec![false; topology.bonds.len()];
    let mut first_positions = BTreeMap::new();
    for (position, bond) in selected.iter().copied().enumerate() {
        if bond.index() >= topology.bonds.len() {
            return Err(MatrixError::ActiveBondOutOfRange {
                position,
                bond,
                bond_count: topology.bonds.len(),
            });
        }
        if let Some(first_position) = first_positions.insert(bond, position) {
            return Err(MatrixError::DuplicateActiveBond {
                bond,
                first_position,
                second_position: position,
            });
        }
        included[bond.index()] = true;
    }
    Ok(Some(included))
}

pub fn adjacency_matrix(
    topology: &TopologyBlock,
    params: &AdjacencyMatrixParams,
) -> Result<DenseMatrix, MatrixError> {
    topology.validate().map_err(MatrixError::InvalidTopology)?;
    let dimension = topology.atoms.len();
    let length = matrix_len(dimension)?;
    let filter = validated_bond_filter(topology, params.bonds_to_use.as_deref())?;

    // BEGIN RDKIT CPP FUNCTION MolOps::getAdjacencyMatrix
    // RDKit✔️✔️: int nAts = mol.getNumAtoms();
    // RDKit✔️✔️: auto *res = new double[nAts * nAts];
    // RDKit✔️✔️: memset(static_cast<void *>(res), emptyVal,
    // RDKit✔️✔️:        nAts * nAts * sizeof(double));
    let empty = f64::from_ne_bytes([params.empty_value as u8; size_of::<f64>()]);
    let mut values = vec![empty; length];

    // RDKit✔️✔️: for (ROMol::ConstBondIterator bondIt = mol.beginBonds();
    // RDKit✔️✔️:      bondIt != mol.endBonds(); bondIt++) {
    for bond in &topology.bonds {
        // RDKit✔️✔️:   if (bondsToUse && !(*bondsToUse)[(*bondIt)->getIdx()]) {
        // RDKit✔️✔️:     continue;
        // RDKit✔️✔️:   }
        if filter
            .as_ref()
            .is_some_and(|included| !included[bond.id().index()])
        {
            continue;
        }
        let begin = bond.begin();
        let end = bond.end();
        if !params.use_bond_order {
            // RDKit✔️✔️:   if (!useBO) {
            // RDKit✔️✔️:     int beg = (*bondIt)->getBeginAtomIdx();
            // RDKit✔️✔️:     int end = (*bondIt)->getEndAtomIdx();
            // RDKit✔️✔️:     res[beg * nAts + end] = 1;
            // RDKit✔️✔️:     res[end * nAts + beg] = 1;
            values[begin.index() * dimension + end.index()] = 1.0;
            values[end.index() * dimension + begin.index()] = 1.0;
        } else {
            // RDKit✔️✔️:   } else {
            // RDKit✔️✔️:     int begIdx = (*bondIt)->getBeginAtomIdx();
            // RDKit✔️✔️:     int endIdx = (*bondIt)->getEndAtomIdx();
            // RDKit✔️✔️:     Atom const *beg = mol.getAtomWithIdx(begIdx);
            // RDKit✔️✔️:     Atom const *end = mol.getAtomWithIdx(endIdx);
            // RDKit✔️✔️:     res[begIdx * nAts + endIdx] =
            // RDKit✔️✔️:         (*bondIt)->getValenceContrib(beg);
            // RDKit✔️✔️:     res[endIdx * nAts + begIdx] =
            // RDKit✔️✔️:         (*bondIt)->getValenceContrib(end);
            values[begin.index() * dimension + end.index()] =
                source_valence_contribution(bond, begin)?;
            values[end.index() * dimension + begin.index()] =
                source_valence_contribution(bond, end)?;
        }
    }
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION MolOps::getAdjacencyMatrix
    Ok(DenseMatrix { dimension, values })
}

fn validate_active_atoms(
    topology: &TopologyBlock,
    active_atoms: &[AtomId],
) -> Result<Vec<Option<usize>>, MatrixError> {
    let mut full_to_active = vec![None; topology.atoms.len()];
    for (position, atom) in active_atoms.iter().copied().enumerate() {
        if atom.index() >= topology.atoms.len() {
            return Err(MatrixError::ActiveAtomOutOfRange {
                position,
                atom,
                atom_count: topology.atoms.len(),
            });
        }
        if let Some(first_position) = full_to_active[atom.index()] {
            return Err(MatrixError::DuplicateActiveAtom {
                atom,
                first_position,
                second_position: position,
            });
        }
        full_to_active[atom.index()] = Some(position);
    }
    Ok(full_to_active)
}

fn selected_distance_bonds<'a>(
    topology: &'a TopologyBlock,
    explicit: Option<&[BondId]>,
    full_to_active: &[Option<usize>],
) -> Result<Vec<&'a Bond>, MatrixError> {
    let filter = validated_bond_filter(topology, explicit)?;
    let mut bonds = Vec::new();
    for bond in &topology.bonds {
        if filter
            .as_ref()
            .is_some_and(|included| !included[bond.id().index()])
        {
            continue;
        }
        let begin_active = full_to_active[bond.begin().index()];
        let end_active = full_to_active[bond.end().index()];
        if filter.is_some() {
            if begin_active.is_none() {
                return Err(MatrixError::ActiveBondEndpointMissing {
                    bond: bond.id(),
                    endpoint: "begin",
                    atom: bond.begin(),
                });
            }
            if end_active.is_none() {
                return Err(MatrixError::ActiveBondEndpointMissing {
                    bond: bond.id(),
                    endpoint: "end",
                    atom: bond.end(),
                });
            }
        }
        if begin_active.is_some() && end_active.is_some() {
            bonds.push(bond);
        }
    }
    Ok(bonds)
}

fn floyd_warshall(
    dimension: usize,
    distances: &mut [f64],
    active_indices: &[usize],
) -> Result<Vec<i32>, MatrixError> {
    let length = matrix_len(dimension)?;
    debug_assert_eq!(distances.len(), length);

    // BEGIN RDKIT CPP FUNCTION FloydWarshall
    // RDKit✔️✔️: currD = new T[dim * dim];
    // RDKit✔️✔️: currP = new int[dim * dim];
    // RDKit✔️✔️: lastD = new T[dim * dim];
    // RDKit✔️✔️: lastP = new int[dim * dim];
    // RDKit✔️✔️: memcpy(static_cast<void *>(lastD), static_cast<void *>(adjMat),
    // RDKit✔️✔️:        dim * dim * sizeof(T));
    let mut last_distances = distances.to_vec();
    let mut current_distances = last_distances.clone();
    let mut paths = vec![0_i32; length];

    // RDKit✔️✔️: for (auto ai : activeAtoms) {
    // RDKit✔️✔️:   int itab = ai * dim;
    // RDKit✔️✔️:   for (int activeAtom : activeAtoms) {
    // RDKit✔️✔️:     if (ai == activeAtom || adjMat[itab + activeAtom] == LOCAL_INF) {
    // RDKit✔️✔️:       pathMat[itab + activeAtom] = -1;
    // RDKit✔️✔️:     } else {
    // RDKit✔️✔️:       pathMat[itab + activeAtom] = ai;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    for &row in active_indices {
        for &column in active_indices {
            paths[row * dimension + column] =
                if row == column || distances[row * dimension + column] == LOCAL_INF {
                    -1
                } else {
                    i32::try_from(row).unwrap_or(i32::MAX)
                };
        }
    }
    let mut last_paths = paths.clone();
    let mut current_paths = last_paths.clone();

    // RDKit✔️✔️: for (auto ak : activeAtoms) {
    // RDKit✔️✔️:   int ktab = ak * dim;
    // RDKit✔️✔️:   for (auto ai : activeAtoms) {
    // RDKit✔️✔️:     int itab = ai * dim;
    // RDKit✔️✔️:     for (int activeAtom : activeAtoms) {
    // RDKit✔️✔️:       T v1 = lastD[itab + activeAtom];
    // RDKit✔️✔️:       T v2 = lastD[itab + ak] + lastD[ktab + activeAtom];
    // RDKit✔️✔️:       if (v1 <= v2) {
    // RDKit✔️✔️:         currD[itab + activeAtom] = v1;
    // RDKit✔️✔️:         currP[itab + activeAtom] = lastP[itab + activeAtom];
    // RDKit✔️✔️:       } else {
    // RDKit✔️✔️:         currD[itab + activeAtom] = v2;
    // RDKit✔️✔️:         currP[itab + activeAtom] = lastP[ktab + activeAtom];
    // RDKit✔️✔️:       }
    for &via in active_indices {
        current_distances.clone_from(&last_distances);
        current_paths.clone_from(&last_paths);
        for &row in active_indices {
            for &column in active_indices {
                let direct = last_distances[row * dimension + column];
                let through = last_distances[row * dimension + via]
                    + last_distances[via * dimension + column];
                if direct <= through {
                    current_distances[row * dimension + column] = direct;
                    current_paths[row * dimension + column] = last_paths[row * dimension + column];
                } else {
                    current_distances[row * dimension + column] = through;
                    current_paths[row * dimension + column] = last_paths[via * dimension + column];
                }
            }
        }
        // RDKit✔️✔️:     }
        // RDKit✔️✔️:   }
        // RDKit✔️✔️:   tTemp = currD;
        // RDKit✔️✔️:   currD = lastD;
        // RDKit✔️✔️:   lastD = tTemp;
        // RDKit✔️✔️:   iTemp = currP;
        // RDKit✔️✔️:   currP = lastP;
        // RDKit✔️✔️:   lastP = iTemp;
        std::mem::swap(&mut current_distances, &mut last_distances);
        std::mem::swap(&mut current_paths, &mut last_paths);
    }
    // RDKit✔️✔️: }
    // RDKit✔️✔️: memcpy(static_cast<void *>(adjMat), static_cast<void *>(lastD),
    // RDKit✔️✔️:        dim * dim * sizeof(T));
    // RDKit✔️✔️: memcpy(static_cast<void *>(pathMat), static_cast<void *>(lastP),
    // RDKit✔️✔️:        dim * dim * sizeof(int));
    // END RDKIT CPP FUNCTION FloydWarshall
    distances.copy_from_slice(&last_distances);
    Ok(last_paths)
}

pub fn topological_distance_matrix(
    topology: &TopologyBlock,
    params: &TopologicalDistanceMatrixParams,
) -> Result<DenseMatrix, MatrixError> {
    topology.validate().map_err(MatrixError::InvalidTopology)?;
    if params.active_atoms.is_none() && params.active_bonds.is_some() {
        return Err(MatrixError::ActiveBondsWithoutAtoms);
    }

    let owned_full_atoms;
    let active_atoms = if let Some(active_atoms) = params.active_atoms.as_deref() {
        active_atoms
    } else {
        owned_full_atoms = (0..topology.atoms.len())
            .map(AtomId::new)
            .collect::<Vec<_>>();
        &owned_full_atoms
    };
    let full_to_active = validate_active_atoms(topology, active_atoms)?;
    let bonds = selected_distance_bonds(topology, params.active_bonds.as_deref(), &full_to_active)?;
    let dimension = active_atoms.len();
    let length = matrix_len(dimension)?;

    // BEGIN RDKIT CPP FUNCTION MolOps::getDistanceMat
    // RDKit✔️✔️: int nAts = mol.getNumAtoms();
    // RDKit✔️✔️: auto *dMat = new double[nAts * nAts];
    // RDKit✔️✔️: for (i = 0; i < nAts * nAts; i++) {
    // RDKit✔️✔️:   dMat[i] = LOCAL_INF;
    // RDKit✔️✔️: }
    // RDKit✔️✔️: for (i = 0; i < nAts; i++) {
    // RDKit✔️✔️:   dMat[i * nAts + i] = 0.0;
    // RDKit✔️✔️: }
    let mut values = vec![LOCAL_INF; length];
    for row in 0..dimension {
        values[row * dimension + row] = 0.0;
    }

    // RDKit✔️✔️: while (firstB != lastB) {
    // RDKit✔️✔️:   const Bond *bond = mol[*firstB];
    // RDKit✔️✔️:   i = bond->getBeginAtomIdx();
    // RDKit✔️✔️:   j = bond->getEndAtomIdx();
    // RDKit✔️✔️:   double contrib;
    for bond in bonds {
        let begin = full_to_active[bond.begin().index()]
            .expect("selected bonds have active begin endpoints");
        let end =
            full_to_active[bond.end().index()].expect("selected bonds have active end endpoints");
        // RDKit✔️✔️:   if (useBO) {
        // RDKit✔️✔️:     if (!bond->getIsAromatic()) {
        // RDKit✔️✔️:       contrib = 1. / bond->getBondTypeAsDouble();
        // RDKit✔️✔️:     } else {
        // RDKit✔️✔️:       contrib = 2. / 3.;
        // RDKit✔️✔️:     }
        // RDKit✔️✔️:   } else {
        // RDKit✔️✔️:     contrib = 1.0;
        // RDKit✔️✔️:   }
        let contribution = if !params.use_bond_order {
            1.0
        } else if bond.is_aromatic() {
            2.0 / 3.0
        } else {
            1.0 / source_bond_order(bond)?
        };
        // RDKit✔️✔️:   dMat[i * nAts + j] = contrib;
        // RDKit✔️✔️:   dMat[j * nAts + i] = contrib;
        // RDKit✔️✔️:   ++firstB;
        // RDKit✔️✔️: }
        values[begin * dimension + end] = contribution;
        values[end * dimension + begin] = contribution;
    }

    let all_rows = (0..dimension).collect::<Vec<_>>();
    let _paths = floyd_warshall(dimension, &mut values, &all_rows)?;

    // RDKit✔️✔️: if (useAtomWts) {
    // RDKit✔️✔️:   for (i = 0; i < nAts; i++) {
    // RDKit✔️✔️:     int anum = mol.getAtomWithIdx(i)->getAtomicNum();
    // RDKit✔️✔️:     dMat[i * nAts + i] = 6.0 / anum;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    if params.use_atom_weights {
        for (row, atom) in active_atoms.iter().copied().enumerate() {
            values[row * dimension + row] =
                6.0 / f64::from(topology.atoms[atom.index()].atomic_number());
        }
    }
    // END RDKIT CPP FUNCTION MolOps::getDistanceMat
    Ok(DenseMatrix { dimension, values })
}

pub(crate) fn unweighted_distance_steps(
    topology: &TopologyBlock,
) -> Result<Vec<usize>, MatrixError> {
    let matrix =
        topological_distance_matrix(topology, &TopologicalDistanceMatrixParams::default())?;
    Ok(matrix
        .values
        .into_iter()
        .map(|distance| {
            if distance == LOCAL_INF {
                usize::MAX
            } else {
                distance as usize
            }
        })
        .collect())
}

pub fn distance_matrix_3d(
    topology: &TopologyBlock,
    coordinates: &CoordinateBlock,
    params: &DistanceMatrix3dParams,
) -> Result<DenseMatrix, MatrixError> {
    topology.validate().map_err(MatrixError::InvalidTopology)?;
    coordinates
        .validate_for_atom_count(topology.atoms.len())
        .map_err(MatrixError::InvalidCoordinates)?;

    // BEGIN RDKIT CPP FUNCTION MolOps::get3DDistanceMat
    // RDKit✔️✔️: const Conformer &conf = mol.getConformer(confId);
    let conformer = match params.conformer_id {
        Some(conformer_id) => coordinates
            .conformers_3d
            .iter()
            .find(|conformer| conformer.id() == conformer_id)
            .ok_or(MatrixError::ConformerNotFound { conformer_id })?,
        None => coordinates
            .conformers_3d
            .first()
            .ok_or(MatrixError::No3dConformer)?,
    };
    // RDKit✔️✔️: unsigned int nAts = mol.getNumAtoms();
    // RDKit✔️✔️: auto *dMat = new double[nAts * nAts];
    let dimension = topology.atoms.len();
    let mut values = vec![0.0; matrix_len(dimension)?];

    // RDKit✔️✔️: for (unsigned int i = 0; i < nAts; ++i) {
    for row in 0..dimension {
        // RDKit✔️✔️:   if (useAtomWts) {
        // RDKit✔️✔️:     dMat[i * nAts + i] =
        // RDKit✔️✔️:         6.0 / mol.getAtomWithIdx(i)->getAtomicNum();
        // RDKit✔️✔️:   } else {
        // RDKit✔️✔️:     dMat[i * nAts + i] = 0.0;
        // RDKit✔️✔️:   }
        values[row * dimension + row] = if params.use_atom_weights {
            6.0 / f64::from(topology.atoms[row].atomic_number())
        } else {
            0.0
        };
        // RDKit✔️✔️:   for (unsigned int j = i + 1; j < nAts; ++j) {
        // RDKit✔️✔️:     double dist =
        // RDKit✔️✔️:         (conf.getAtomPos(i) - conf.getAtomPos(j)).length();
        // RDKit✔️✔️:     dMat[i * nAts + j] = dist;
        // RDKit✔️✔️:     dMat[j * nAts + i] = dist;
        for column in row + 1..dimension {
            let left = conformer.coordinates()[row];
            let right = conformer.coordinates()[column];
            let dx = left[0] - right[0];
            let dy = left[1] - right[1];
            let dz = left[2] - right[2];
            let distance = (dx * dx + dy * dy + dz * dz).sqrt();
            values[row * dimension + column] = distance;
            values[column * dimension + row] = distance;
        }
        // RDKit✔️✔️:   }
    }
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION MolOps::get3DDistanceMat
    Ok(DenseMatrix { dimension, values })
}

#[cfg(test)]
mod tests {
    use super::{LOCAL_INF, floyd_warshall};

    #[test]
    fn floyd_warshall_preserves_ties_and_limits_updates_to_ordered_active_indices() {
        let dimension = 5;
        let mut distances = vec![77.0; dimension * dimension];
        for &row in &[1_usize, 2, 3, 4] {
            for &column in &[1_usize, 2, 3, 4] {
                distances[row * dimension + column] = if row == column { 0.0 } else { LOCAL_INF };
            }
        }
        for (left, right) in [(1_usize, 2_usize), (1, 3), (2, 4), (3, 4)] {
            distances[left * dimension + right] = 1.0;
            distances[right * dimension + left] = 1.0;
        }

        let paths = floyd_warshall(dimension, &mut distances, &[1, 2, 3, 4]).unwrap();
        assert_eq!(distances[1 * dimension + 4], 2.0);
        assert_eq!(distances[4 * dimension + 1], 2.0);
        assert_eq!(paths[1 * dimension + 4], 2);
        assert_eq!(paths[4 * dimension + 1], 2);
        assert_eq!(&distances[..dimension], &[77.0; 5]);

        let mut active_only = vec![91.0; dimension * dimension];
        for &row in &[1_usize, 3, 4] {
            for &column in &[1_usize, 3, 4] {
                active_only[row * dimension + column] = if row == column { 0.0 } else { LOCAL_INF };
            }
        }
        for (left, right) in [(1_usize, 3_usize), (3, 4)] {
            active_only[left * dimension + right] = 1.0;
            active_only[right * dimension + left] = 1.0;
        }
        floyd_warshall(dimension, &mut active_only, &[1, 3, 4]).unwrap();
        assert_eq!(active_only[1 * dimension + 4], 2.0);
        assert_eq!(active_only[4 * dimension + 1], 2.0);
        assert_eq!(&active_only[..dimension], &[91.0; 5]);
    }
}
