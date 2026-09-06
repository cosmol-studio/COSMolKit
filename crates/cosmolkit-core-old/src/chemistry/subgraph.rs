//! Shared source-reproduced molecular subgraph traversal.
//!
//! Fingerprints and graph descriptors consume this single RDKit path-ordering
//! implementation. This module does not own descriptor or fingerprint policy.

use std::collections::BTreeMap;

use crate::Molecule;

#[cfg(test)]
mod tests;

#[derive(Debug, Clone, PartialEq, Eq, thiserror::Error)]
pub(crate) enum SubgraphPathError {
    #[error("invalid subgraph path arguments: {reason}")]
    InvalidArguments { reason: &'static str },
}

pub(crate) fn extend_paths(
    adjacency: &[u8],
    dim: usize,
    paths: &[Vec<usize>],
    allow_ring_closures: i64,
    distance_matrix: Option<&[f64]>,
) -> Result<Vec<Vec<usize>>, SubgraphPathError> {
    // BEGIN RDKIT CPP FUNCTION Subgraphs::extendPaths
    // RDKit✔️✔️: PATH_LIST
    // RDKit✔️✔️: extendPaths(int *adjMat, unsigned int dim, const PATH_LIST &paths,
    // RDKit✔️✔️:             int allowRingClosures = -1, double *distMat = nullptr) {
    // RDKit✔️✔️:   PRECONDITION(adjMat, "no matrix");
    // RDKit✔️✔️:   //
    // RDKit✔️✔️:   //  extend each of the currently active paths by adding
    // RDKit✔️✔️:   //   a single adjacent index to the end of each
    // RDKit✔️✔️:   //
    // RDKit✔️✔️:   PATH_LIST res;
    // RDKit✔️✔️:   PATH_LIST::const_iterator path;
    // RDKit✔️✔️:   for (path = paths.begin(); path != paths.end(); ++path) {
    // RDKit✔️✔️:     unsigned int endIdx = (*path)[path->size() - 1];
    // RDKit✔️✔️:     unsigned int iTab = endIdx * dim;
    // RDKit✔️✔️:     for (unsigned int otherIdx = 0; otherIdx < dim; otherIdx++) {
    // RDKit✔️✔️:       if (adjMat[iTab + otherIdx] == 1) {
    // RDKit✔️✔️:         if (distMat &&
    // RDKit✔️✔️:             distMat[path->front() * dim + otherIdx] - path->size() < -0.001) {
    // RDKit✔️✔️:           continue;
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:         // test 1: make sure the new atom is not already
    // RDKit✔️✔️:         //   in the path
    // RDKit✔️✔️:         auto loc =
    // RDKit✔️✔️:             std::find(path->begin(), path->end(), static_cast<int>(otherIdx));
    // RDKit✔️✔️:         // The two conditions for adding the atom are:
    // RDKit✔️✔️:         //   1) it's not there already
    // RDKit✔️✔️:         //   2) it's there, but ring closures are allowed and this
    // RDKit✔️✔️:         //      will be the last addition to the path.
    // RDKit✔️✔️:         if (loc == path->end()) {
    // RDKit✔️✔️:           // the easy case
    // RDKit✔️✔️:           // PATH_TYPE newPath=*path;
    // RDKit✔️✔️:           // newPath.push_back(otherIdx);
    // RDKit✔️✔️:           // res.push_back(newPath);
    // RDKit✔️✔️:           res.push_back(*path);
    // RDKit✔️✔️:           res.rbegin()->push_back(otherIdx);
    // RDKit✔️✔️:         } else if (allowRingClosures > 2 &&
    // RDKit✔️✔️:                    static_cast<int>(path->size()) == allowRingClosures - 1) {
    // RDKit✔️✔️:           // We *might* be adding the atom, but we need to make sure
    // RDKit✔️✔️:           // that we're not just duplicating the second to last
    // RDKit✔️✔️:           // element of the path:
    // RDKit✔️✔️:           auto rIt = path->rbegin();
    // RDKit✔️✔️:           rIt++;
    // RDKit✔️✔️:           if (*rIt != static_cast<int>(otherIdx)) {
    // RDKit✔️✔️:             // PATH_TYPE newPath=*path;
    // RDKit✔️✔️:             // newPath.push_back(otherIdx);
    // RDKit✔️✔️:             // res.push_back(newPath);
    // RDKit✔️✔️:             res.push_back(*path);
    // RDKit✔️✔️:             res.rbegin()->push_back(otherIdx);
    // RDKit✔️✔️:           }
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return res;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION Subgraphs::extendPaths
    let matrix_len = dim
        .checked_mul(dim)
        .ok_or(SubgraphPathError::InvalidArguments {
            reason: "path adjacency matrix dimensions overflow",
        })?;
    if adjacency.len() != matrix_len {
        return Err(SubgraphPathError::InvalidArguments {
            reason: "path adjacency matrix has invalid dimensions",
        });
    }
    if distance_matrix.is_some_and(|matrix| matrix.len() != matrix_len) {
        return Err(SubgraphPathError::InvalidArguments {
            reason: "path distance matrix has invalid dimensions",
        });
    }

    let mut result = Vec::new();
    for path in paths {
        let end = *path.last().ok_or(SubgraphPathError::InvalidArguments {
            reason: "path to extend must not be empty",
        })?;
        let start = path[0];
        if end >= dim || start >= dim {
            return Err(SubgraphPathError::InvalidArguments {
                reason: "path atom index is out of range",
            });
        }
        let row_offset = end * dim;
        // The source scans every column of the dense adjacency matrix in
        // atom-index order. Keep that order rather than the molecule's bond
        // insertion order.
        for other_idx in 0..dim {
            if adjacency[row_offset + other_idx] != 1 {
                continue;
            }
            if distance_matrix.is_some_and(|matrix| {
                matrix[start * dim + other_idx] - (path.len() as f64) < -0.001
            }) {
                continue;
            }
            if !path.contains(&other_idx) {
                let mut next = path.clone();
                next.push(other_idx);
                result.push(next);
            } else if allow_ring_closures > 2
                && i64::try_from(path.len()).ok() == Some(allow_ring_closures - 1)
                && path[path.len() - 2] != other_idx
            {
                let mut next = path.clone();
                next.push(other_idx);
                result.push(next);
            }
        }
    }
    Ok(result)
}

fn path_finder_helper(
    adjacency: &[u8],
    dim: usize,
    min_len: usize,
    max_len: usize,
    rooted_at_atom: i64,
    distance_matrix: Option<&[f64]>,
) -> Result<BTreeMap<usize, Vec<Vec<usize>>>, SubgraphPathError> {
    // BEGIN RDKIT CPP FUNCTION Subgraphs::pathFinderHelper
    // RDKit✔️✔️: INT_PATH_LIST_MAP
    // RDKit✔️✔️: pathFinderHelper(int *adjMat, unsigned int dim, unsigned int minLen,
    // RDKit✔️✔️:                  unsigned int maxLen, int rootedAtAtom, double *distMat) {
    // RDKit✔️✔️:   PRECONDITION(adjMat, "no matrix");
    // RDKit✔️✔️:   PRECONDITION(minLen <= maxLen, "bad lengths provided");
    // RDKit✔️✔️:   // finds all paths of length N using an adjacency matrix,
    // RDKit✔️✔️:   //  which is constructed elsewhere
    // RDKit✔️✔️:   INT_PATH_LIST_MAP res;
    // RDKit✔️✔️:   PATH_LIST paths;
    // RDKit✔️✔️:   paths.clear();
    // RDKit✔️✔️:
    // RDKit✔️✔️:   if (rootedAtAtom < 0) {
    // RDKit✔️✔️:     // start a path at each possible index
    // RDKit✔️✔️:     for (unsigned int i = 0; i < dim; i++) {
    // RDKit✔️✔️:       PATH_TYPE tPath;
    // RDKit✔️✔️:       tPath.push_back(i);
    // RDKit✔️✔️:       paths.push_back(tPath);
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   } else if (rootedAtAtom < static_cast<int>(dim)) {
    // RDKit✔️✔️:     // only start a path at the atom of interest:
    // RDKit✔️✔️:     PATH_TYPE tPath;
    // RDKit✔️✔️:     tPath.push_back(rootedAtAtom);
    // RDKit✔️✔️:     paths.push_back(tPath);
    // RDKit✔️✔️:   } else {
    // RDKit✔️✔️:     return res;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   // and build them up one index at a time:
    // RDKit✔️✔️:   for (unsigned int length = 1; length < maxLen; length++) {
    // RDKit✔️✔️:     // extend each path:
    // RDKit✔️✔️:     if (length >= minLen) {
    // RDKit✔️✔️:       res[length] = paths;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     paths = extendPaths(adjMat, dim, paths, maxLen, distMat);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   res[maxLen] = paths;
    // RDKit✔️✔️:
    // RDKit✔️✔️:   return res;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION Subgraphs::pathFinderHelper
    if min_len > max_len {
        return Err(SubgraphPathError::InvalidArguments {
            reason: "minimum path length exceeds maximum path length",
        });
    }
    let allow_ring_closures =
        i64::try_from(max_len).map_err(|_| SubgraphPathError::InvalidArguments {
            reason: "path length exceeds supported range",
        })?;

    let mut result = BTreeMap::new();
    let mut paths = Vec::new();
    if rooted_at_atom < 0 {
        for atom_index in 0..dim {
            paths.push(vec![atom_index]);
        }
    } else if usize::try_from(rooted_at_atom).is_ok_and(|root| root < dim) {
        paths.push(vec![rooted_at_atom as usize]);
    } else {
        return Ok(result);
    }

    for length in 1..max_len {
        if length >= min_len {
            result.insert(length, paths.clone());
        }
        paths = extend_paths(adjacency, dim, &paths, allow_ring_closures, distance_matrix)?;
    }
    result.insert(max_len, paths);
    Ok(result)
}

pub(crate) fn rdkit_fp_bond_between_atoms(
    molecule: &Molecule,
    begin: usize,
    end: usize,
) -> Option<usize> {
    // RDKit source: GraphMol/ROMol.cpp lines 338-350.
    // RDKit✔️✔️: const Bond *ROMol::getBondBetweenAtoms(unsigned int idx1,
    // RDKit✔️✔️:                                        unsigned int idx2) const {
    // RDKit✔️✔️:   URANGE_CHECK(idx1, getNumAtoms());
    // RDKit✔️✔️:   URANGE_CHECK(idx2, getNumAtoms());
    // RDKit✔️✔️:   const Bond *res = nullptr;
    // RDKit✔️✔️:   auto [edge, found] = boost::edge(boost::vertex(idx1, d_graph),
    // RDKit✔️✔️:                                    boost::vertex(idx2, d_graph), d_graph);
    // RDKit✔️✔️:   if (found) {
    // RDKit✔️✔️:     res = d_graph[edge];
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return res;
    // RDKit✔️✔️: }
    if begin >= molecule.num_atoms() || end >= molecule.num_atoms() {
        return None;
    }
    molecule
        .topology_block()
        .adjacency
        .neighbors_of(begin)
        .iter()
        .find(|neighbor| neighbor.atom_index == end)
        .map(|neighbor| neighbor.bond.index())
}

pub(crate) fn find_all_paths_of_lengths_m_to_n(
    molecule: &Molecule,
    mut lower_len: usize,
    mut upper_len: usize,
    use_bonds: bool,
    use_hs: bool,
    rooted_at_atom: i64,
    only_shortest_paths: bool,
) -> Result<BTreeMap<usize, Vec<Vec<usize>>>, SubgraphPathError> {
    // BEGIN RDKIT CPP FUNCTION findAllPathsOfLengthsMtoN
    // RDKit✔️✔️: INT_PATH_LIST_MAP
    // RDKit✔️✔️: findAllPathsOfLengthsMtoN(const ROMol &mol, unsigned int lowerLen,
    // RDKit✔️✔️:                           unsigned int upperLen, bool useBonds, bool useHs,
    // RDKit✔️✔️:                           int rootedAtAtom, bool onlyShortestPaths) {
    // RDKit✔️✔️:   //
    // RDKit✔️✔️:   //  We can't be clever here and just use the bond adjacency matrix
    // RDKit✔️✔️:   //  to solve this problem when useBonds is true.  This is because
    // RDKit✔️✔️:   //  the bond adjacency matrices for the molecules C1CC1 and CC(C)C
    // RDKit✔️✔️:   //  are indistinguishable.  In the second case, t-butane (and
    // RDKit✔️✔️:   //  anything else with a T junction), we'll get some subgraphs mixed
    // RDKit✔️✔️:   //  in with the paths.  So we have to construct paths of atoms and
    // RDKit✔️✔️:   //  then convert them into bond paths.
    // RDKit✔️✔️:   //
    // RDKit✔️✔️:   PRECONDITION(lowerLen <= upperLen, "");
    // RDKit✔️✔️:
    // RDKit✔️✔️:   // the molecule owns the distance matrix pointer (if we need to get it)
    // RDKit✔️✔️:   double *distMat = onlyShortestPaths ? MolOps::getDistanceMat(mol) : nullptr;
    // RDKit✔️✔️:   int *adjMat, dim;
    // RDKit✔️✔️:   dim = mol.getNumAtoms();
    // RDKit✔️✔️:   adjMat = new int[dim * dim];
    // RDKit✔️✔️:   memset((void *)adjMat, 0, dim * dim * sizeof(int));
    // RDKit✔️✔️:
    // RDKit✔️✔️:   if (!distMat) {
    // RDKit✔️✔️:     // generate the adjacency matrix by hand by looping over the bonds
    // RDKit✔️✔️:     ROMol::ConstBondIterator bondIt;
    // RDKit✔️✔️:     for (bondIt = mol.beginBonds(); bondIt != mol.endBonds(); bondIt++) {
    // RDKit✔️✔️:       Atom *beg = (*bondIt)->getBeginAtom();
    // RDKit✔️✔️:       Atom *end = (*bondIt)->getEndAtom();
    // RDKit✔️✔️:       // check for H, which we might be skipping
    // RDKit✔️✔️:       if (useHs || (beg->getAtomicNum() != 1 && end->getAtomicNum() != 1)) {
    // RDKit✔️✔️:         adjMat[beg->getIdx() * dim + end->getIdx()] = 1;
    // RDKit✔️✔️:         adjMat[end->getIdx() * dim + beg->getIdx()] = 1;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   } else {
    // RDKit✔️✔️:     // if we have the distance matrix, we can just loop over that:
    // RDKit✔️✔️:     for (auto i = 0; i < dim; ++i) {
    // RDKit✔️✔️:       for (auto j = i + 1; j < dim; ++j) {
    // RDKit✔️✔️:         if (fabs(distMat[i * dim + j] - 1) < 1e-4) {
    // RDKit✔️✔️:           adjMat[i * dim + j] = 1;
    // RDKit✔️✔️:           adjMat[j * dim + i] = 1;
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   // if we're using bonds, we'll need to find paths of length N+1,
    // RDKit✔️✔️:   // then convert them
    // RDKit✔️✔️:   if (useBonds) {
    // RDKit✔️✔️:     ++lowerLen;
    // RDKit✔️✔️:     ++upperLen;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   // find the paths themselves
    // RDKit✔️✔️:   INT_PATH_LIST_MAP atomPaths = Subgraphs::pathFinderHelper(
    // RDKit✔️✔️:       adjMat, dim, lowerLen, upperLen, rootedAtAtom, distMat);
    // RDKit✔️✔️:
    // RDKit✔️✔️:   // clean up the adjacency matrix
    // RDKit✔️✔️:   delete[] adjMat;
    // RDKit✔️✔️:
    // RDKit✔️✔️:   INT_PATH_LIST_MAP res;
    // RDKit✔️✔️:
    // RDKit✔️✔️:   //
    // RDKit✔️✔️:   //--------------------------------------------------------
    // RDKit✔️✔️:   // loop through all the paths we have and make sure that there are
    // RDKit✔️✔️:   // no duplicates (duplicate = contains identical bond indices)
    // RDKit✔️✔️:   //
    // RDKit✔️✔️:   //  We need to use the bond paths for this duplicate finding
    // RDKit✔️✔️:   //  because, in rings, there can be many paths which share atom
    // RDKit✔️✔️:   //  indices but which have different bond compositions. For example,
    // RDKit✔️✔️:   //  there is only one "atom unique" path of length 5 bonds (6 atoms)
    // RDKit✔️✔️:   //  through a 6-ring, but there are six bond paths.
    // RDKit✔️✔️:   //
    // RDKit✔️✔️:   if (!useBonds && lowerLen >= 1) {
    // RDKit✔️✔️:     res[1] = atomPaths[1];
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   if (useBonds || upperLen > 1) {
    // RDKit✔️✔️:     for (unsigned int i = lowerLen; i <= upperLen; ++i) {
    // RDKit✔️✔️:       if (i <= 1) {
    // RDKit✔️✔️:         continue;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:
    // RDKit✔️✔️:       std::vector<boost::dynamic_bitset<>> invars;
    // RDKit✔️✔️:
    // RDKit✔️✔️:       for (PATH_LIST::const_iterator vivI = atomPaths[i].begin();
    // RDKit✔️✔️:            vivI != atomPaths[i].end(); ++vivI) {
    // RDKit✔️✔️:         boost::dynamic_bitset<> invar(mol.getNumBonds());
    // RDKit✔️✔️:         const PATH_TYPE &resi = *vivI;
    // RDKit✔️✔️:         PATH_TYPE locV;
    // RDKit✔️✔️:         locV.reserve(i);
    // RDKit✔️✔️:         for (unsigned int j = 0; j < i - 1; j++) {
    // RDKit✔️✔️:           const Bond *bond = mol.getBondBetweenAtoms(resi[j], resi[j + 1]);
    // RDKit✔️✔️:           locV.push_back(bond->getIdx());
    // RDKit✔️✔️:           invar.set(bond->getIdx());
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:         if (std::find(invars.begin(), invars.end(), invar) == invars.end()) {
    // RDKit✔️✔️:           invars.push_back(invar);
    // RDKit✔️✔️:           if (useBonds) {
    // RDKit✔️✔️:             res[i - 1].push_back(locV);
    // RDKit✔️✔️:           } else {
    // RDKit✔️✔️:             res[i].push_back(resi);
    // RDKit✔️✔️:           }
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return res;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION findAllPathsOfLengthsMtoN
    if lower_len > upper_len {
        return Err(SubgraphPathError::InvalidArguments {
            reason: "minimum path length exceeds maximum path length",
        });
    }

    let distance_matrix = only_shortest_paths
        .then(|| crate::chemistry::matrices::topological_distance_matrix(molecule));
    let dim = molecule.num_atoms();
    let matrix_len = dim
        .checked_mul(dim)
        .ok_or(SubgraphPathError::InvalidArguments {
            reason: "path adjacency matrix dimensions overflow",
        })?;
    let mut adjacency = vec![0u8; matrix_len];
    if let Some(matrix) = distance_matrix.as_deref() {
        for first in 0..dim {
            for second in (first + 1)..dim {
                if (matrix[first * dim + second] - 1.0).abs() < 1.0e-4 {
                    adjacency[first * dim + second] = 1;
                    adjacency[second * dim + first] = 1;
                }
            }
        }
    } else {
        for bond in molecule.bonds() {
            let begin = bond.begin().index();
            let end = bond.end().index();
            if use_hs
                || (molecule.atoms()[begin].atomic_number() != 1
                    && molecule.atoms()[end].atomic_number() != 1)
            {
                adjacency[begin * dim + end] = 1;
                adjacency[end * dim + begin] = 1;
            }
        }
    }

    if use_bonds {
        lower_len = lower_len
            .checked_add(1)
            .ok_or(SubgraphPathError::InvalidArguments {
                reason: "minimum bond path length exceeds supported range",
            })?;
        upper_len = upper_len
            .checked_add(1)
            .ok_or(SubgraphPathError::InvalidArguments {
                reason: "maximum bond path length exceeds supported range",
            })?;
    }
    let atom_paths = path_finder_helper(
        &adjacency,
        dim,
        lower_len,
        upper_len,
        rooted_at_atom,
        distance_matrix.as_deref(),
    )?;

    let mut result = BTreeMap::new();
    if !use_bonds && lower_len >= 1 {
        result.insert(1, atom_paths.get(&1).cloned().unwrap_or_default());
    }
    if use_bonds || upper_len > 1 {
        for path_length in lower_len..=upper_len {
            if path_length <= 1 {
                continue;
            }
            let mut invariants: Vec<Vec<bool>> = Vec::new();
            if let Some(paths) = atom_paths.get(&path_length) {
                for atom_path in paths {
                    if atom_path.len() < path_length {
                        return Err(SubgraphPathError::InvalidArguments {
                            reason: "enumerated atom path is shorter than its length key",
                        });
                    }
                    let mut invariant = vec![false; molecule.num_bonds()];
                    let mut bond_path = Vec::with_capacity(path_length);
                    for atom_pair in atom_path[..path_length].windows(2) {
                        let bond_index =
                            rdkit_fp_bond_between_atoms(molecule, atom_pair[0], atom_pair[1])
                                .ok_or(SubgraphPathError::InvalidArguments {
                                    reason: "path contains no connecting bond",
                                })?;
                        bond_path.push(bond_index);
                        invariant[bond_index] = true;
                    }
                    if !invariants.contains(&invariant) {
                        invariants.push(invariant);
                        if use_bonds {
                            result.entry(path_length - 1).or_default().push(bond_path);
                        } else {
                            result
                                .entry(path_length)
                                .or_default()
                                .push(atom_path.clone());
                        }
                    }
                }
            }
        }
    }
    Ok(result)
}

pub(crate) fn find_all_paths_of_length_n(
    molecule: &Molecule,
    target_len: usize,
    use_bonds: bool,
    use_hs: bool,
    rooted_at_atom: i64,
    only_shortest_paths: bool,
) -> Result<Vec<Vec<usize>>, SubgraphPathError> {
    // BEGIN RDKIT CPP FUNCTION findAllPathsOfLengthN
    // RDKit✔️✔️: PATH_LIST
    // RDKit✔️✔️: findAllPathsOfLengthN(const ROMol &mol, unsigned int targetLen, bool useBonds,
    // RDKit✔️✔️:                       bool useHs, int rootedAtAtom, bool onlyShortestPaths) {
    // RDKit✔️✔️:   return findAllPathsOfLengthsMtoN(mol, targetLen, targetLen, useBonds, useHs,
    // RDKit✔️✔️:                                    rootedAtAtom, onlyShortestPaths)[targetLen];
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION findAllPathsOfLengthN
    let mut paths = find_all_paths_of_lengths_m_to_n(
        molecule,
        target_len,
        target_len,
        use_bonds,
        use_hs,
        rooted_at_atom,
        only_shortest_paths,
    )?;
    Ok(paths.remove(&target_len).unwrap_or_default())
}
