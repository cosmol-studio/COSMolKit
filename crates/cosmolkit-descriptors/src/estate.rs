//! Electrotopological state of the atoms in a supplied molecular graph.
//! RDKit 4ce8a0636f9b8ef28d6c3f4af4f502392c997d6b, GraphMol/Descriptors/EState.cpp.
//! Copyright (C) 2026 Glysade LLC and other RDKit contributors (BSD).

use cosmolkit_core::{
    PathError, ValenceAssignment, element_info, topological_distance_matrix,
    total_hydrogen_count_from_cache,
};
use cosmolkit_model::TopologyBlock;

use crate::{DescriptorError, DescriptorResult};

fn principal_quantum_number(atomic_number: u8) -> f64 {
    // RDKit❗✔️: unsigned int principalQuantumNumber(unsigned int atNum) {
    // RDKit❗✔️:   if (atNum <= 2) {
    // RDKit❗✔️:     return 1;
    // RDKit❗✔️:   } else if (atNum <= 10) {
    // RDKit❗✔️:     return 2;
    // RDKit❗✔️:   } else if (atNum <= 18) {
    // RDKit❗✔️:     return 3;
    // RDKit❗✔️:   } else if (atNum <= 36) {
    // RDKit❗✔️:     return 4;
    // RDKit❗✔️:   } else if (atNum <= 54) {
    // RDKit❗✔️:     return 5;
    // RDKit❗✔️:   } else if (atNum <= 86) {
    // RDKit❗✔️:     return 6;
    // RDKit❗✔️:   }
    // RDKit❗✔️:   return 7;
    // RDKit❗✔️: }
    match atomic_number {
        0..=2 => 1.0,
        3..=10 => 2.0,
        11..=18 => 3.0,
        19..=36 => 4.0,
        37..=54 => 5.0,
        55..=86 => 6.0,
        _ => 7.0,
    }
}

/// Atom-ordered indices; only atoms of positive degree read hydrogen counts.
pub fn estate_indices(
    topology: &TopologyBlock,
    valence: Option<&ValenceAssignment>,
) -> DescriptorResult<Vec<f64>> {
    // RDKit❗✔️: std::vector<double> calcEStateIndices(const ROMol &mol) {
    // RDKit❗✔️:   const unsigned int nAtoms = mol.getNumAtoms();
    // RDKit❗✔️:   std::vector<double> Is(nAtoms, 0.0);
    // RDKit❗✔️:   if (!nAtoms) {
    // RDKit❗✔️:     return Is;
    // RDKit❗✔️:   }
    // RDKit❗✔️:   const auto *tbl = PeriodicTable::getTable();
    // RDKit❗✔️:   for (const auto at : mol.atoms()) {
    // RDKit❗✔️:     auto d = at->getDegree();
    // RDKit❗✔️:     if (d > 0) {
    // RDKit❗✔️:       auto atNum = at->getAtomicNum();
    // RDKit❗✔️:       double dv = static_cast<int>(tbl->getNouterElecs(atNum)) -
    // RDKit❗✔️:                   static_cast<int>(at->getTotalNumHs());
    // RDKit❗✔️:       double N = principalQuantumNumber(atNum);
    // RDKit❗✔️:       Is[at->getIdx()] = (4. / (N * N) * dv + 1) / d;
    // RDKit❗✔️:     }
    // RDKit❗✔️:   }
    let n = topology.atoms.len();
    if n == 0 {
        return Ok(Vec::new());
    }

    let mut intrinsic = vec![0.0; n];
    for (i, atom) in topology.atoms.iter().enumerate() {
        let degree = topology.adjacency.neighbors_of(i).len();
        if degree == 0 {
            continue;
        }
        let hydrogens = total_hydrogen_count_from_cache(topology, valence, atom.id(), false)
            .map_err(|source| DescriptorError::Valence {
                function: "estate_indices",
                source,
            })?;
        let dv = f64::from(element_info(atom.element()).outer_electrons) - f64::from(hydrogens);
        let shell = principal_quantum_number(atom.atomic_number());
        intrinsic[i] = (4.0 / (shell * shell) * dv + 1.0) / degree as f64;
    }

    // RDKit❗✔️:   const bool useBO = false;
    // RDKit❗✔️:   const bool useAtomWts = false;
    // RDKit❗✔️:   const bool force = true;
    // RDKit❗✔️:   const double *dists = MolOps::getDistanceMat(mol, useBO, useAtomWts, force);
    // RDKit❗✔️:   std::vector<double> accum(nAtoms, 0.0);
    // RDKit❗✔️:   for (unsigned int i = 0; i < nAtoms; ++i) {
    // RDKit❗✔️:     for (unsigned int j = i + 1; j < nAtoms; ++j) {
    // RDKit❗✔️:       double p = dists[i * nAtoms + j] + 1;
    // RDKit❗✔️:       if (p < 1e6) {
    // RDKit❗✔️:         double tmp = (Is[i] - Is[j]) / (p * p);
    // RDKit❗✔️:         accum[i] += tmp;
    // RDKit❗✔️:         accum[j] -= tmp;
    // RDKit❗✔️:       }
    // RDKit❗✔️:     }
    // RDKit❗✔️:   }
    let distances =
        topological_distance_matrix(topology, &Default::default()).map_err(|source| {
            DescriptorError::Path {
                function: "estate_indices",
                source: PathError::Matrix(source),
            }
        })?;
    let mut indices = vec![0.0; n];
    for i in 0..n {
        for j in i + 1..n {
            let p = distances.values()[i * n + j] + 1.0;
            if p < 1e6 {
                let transfer = (intrinsic[i] - intrinsic[j]) / (p * p);
                indices[i] += transfer;
                indices[j] -= transfer;
            }
        }
    }

    // RDKit❗✔️:   for (unsigned int i = 0; i < nAtoms; ++i) {
    // RDKit❗✔️:     accum[i] += Is[i];
    // RDKit❗✔️:   }
    // RDKit❗✔️:   return accum;
    // RDKit❗✔️: }
    for (index, base) in indices.iter_mut().zip(intrinsic) {
        *index += base;
    }
    Ok(indices)
}

/// Extrema of the signed and absolute values of the same atom-ordered indices.
#[derive(Clone, Copy, Default, Debug, PartialEq)]
pub struct EStateExtrema {
    pub maximum: f64,
    pub minimum: f64,
    pub maximum_absolute: f64,
    pub minimum_absolute: f64,
}

impl EStateExtrema {
    pub fn from_indices(indices: &[f64]) -> Self {
        // RDKit❗✔️: EStates calcEStates(const ROMol &mol) {
        // RDKit❗✔️:   auto estates = calcEStateIndices(mol);
        // RDKit❗✔️:
        // RDKit❗✔️:   if(estates.empty()) {
        // RDKit❗✔️:     return {};
        // RDKit❗✔️:   }
        // RDKit❗✔️:
        // RDKit❗✔️:   EStates res { estates, estates[0], estates[0], std::fabs(estates[0]), std::fabs(estates[0]) };
        // RDKit❗✔️:
        // RDKit❗✔️:   for(auto estate : estates) {
        // RDKit❗✔️:     res.maxEState = std::max(res.maxEState, estate);
        // RDKit❗✔️:     res.minEState = std::min(res.minEState, estate);
        // RDKit❗✔️:     estate = std::fabs(estate);
        // RDKit❗✔️:     res.maxAbsEState = std::max(res.maxAbsEState, estate);
        // RDKit❗✔️:     res.minAbsEState = std::min(res.minAbsEState, estate);
        // RDKit❗✔️:   }
        // RDKit❗✔️:   return res;
        // RDKit❗✔️: }
        let Some((&first, rest)) = indices.split_first() else {
            return Self::default();
        };
        let mut extrema = Self {
            maximum: first,
            minimum: first,
            maximum_absolute: first.abs(),
            minimum_absolute: first.abs(),
        };
        for &index in rest {
            if index > extrema.maximum {
                extrema.maximum = index;
            }
            if index < extrema.minimum {
                extrema.minimum = index;
            }
            let magnitude = index.abs();
            if magnitude > extrema.maximum_absolute {
                extrema.maximum_absolute = magnitude;
            }
            if magnitude < extrema.minimum_absolute {
                extrema.minimum_absolute = magnitude;
            }
        }
        extrema
    }
}

/// Direction of the projection: the first property selects the bin; the second is summed.
#[derive(Clone, Copy, Debug)]
pub enum EStateVsaFamily {
    AreaByEState,
    EStateByArea,
}

impl EStateVsaFamily {
    /// Default vectors refresh Labute; numbered projections read through their own family slot.
    pub fn bins<'a>(
        self,
        topology: &TopologyBlock,
        valence: Option<&ValenceAssignment>,
        state: &'a mut crate::DescriptorComputedState,
        refresh: bool,
    ) -> DescriptorResult<&'a [f64]> {
        // RDKit❗✔️: const std::vector<double> estateBins = {-0.390, 0.290, 0.717, 1.165, 1.540,
        // RDKit❗✔️:                                         1.807,  2.05,  4.69,  9.17,  15.0};
        // RDKit❗✔️: const std::vector<double> vsaEStateBins = {4.78, 5.00, 5.410, 5.740, 6.00,
        // RDKit❗✔️:                                            6.07, 6.45, 7.00,  11.0};
        const ESTATE_CUTOFFS: &[f64] = &[
            -0.390, 0.290, 0.717, 1.165, 1.540, 1.807, 2.05, 4.69, 9.17, 15.0,
        ];
        const AREA_CUTOFFS: &[f64] = &[4.78, 5.00, 5.410, 5.740, 6.00, 6.07, 6.45, 7.00, 11.0];

        if refresh || self.slot(state).is_none() {
            let indices = estate_indices(topology, valence)?;
            // RDKit❗✔️: std::vector<double> labuteContribs(const ROMol &mol) {
            // RDKit❗✔️:   std::vector<double> Vi(mol.getNumAtoms(), 0.0);
            // RDKit❗✔️:   double hContrib = 0.0;
            // RDKit❗✔️:   getLabuteAtomContribs(mol, Vi, hContrib, true, true);
            // RDKit❗✔️:   return Vi;
            // RDKit❗✔️: }
            let areas = crate::labute_atom_rows(topology, state, true, true)?.atom_contributions;
            let result = match self {
                Self::AreaByEState => {
                    // RDKit❗✔️: std::vector<double> calcEState_VSA(const ROMol &mol,
                    // RDKit❗✔️:                                    const std::vector<double> *bins) {
                    // RDKit❗✔️:   if (!bins) {
                    // RDKit❗✔️:     bins = &estateBins;
                    // RDKit❗✔️:   }
                    // RDKit❗✔️:   auto propContribs = calcEStateIndices(mol);
                    // RDKit❗✔️:   auto volContribs = labuteContribs(mol);
                    // RDKit❗✔️:   std::vector<double> res(bins->size() + 1, 0.0);
                    // RDKit❗✔️:   for (unsigned int i = 0; i < propContribs.size(); ++i) {
                    // RDKit❗✔️:     auto nbin = std::upper_bound(bins->begin(), bins->end(), propContribs[i]) -
                    // RDKit❗✔️:                 bins->begin();
                    // RDKit❗✔️:     res[nbin] += volContribs[i];
                    // RDKit❗✔️:   }
                    // RDKit❗✔️:   return res;
                    // RDKit❗✔️: }
                    crate::assign_contribs_to_bins(&areas, &indices, ESTATE_CUTOFFS)?
                }
                Self::EStateByArea => {
                    // RDKit❗✔️: std::vector<double> calcVSA_EState(const ROMol &mol,
                    // RDKit❗✔️:                                    const std::vector<double> *bins) {
                    // RDKit❗✔️:   if (!bins) {
                    // RDKit❗✔️:     bins = &vsaEStateBins;
                    // RDKit❗✔️:   }
                    // RDKit❗✔️:   auto propContribs = calcEStateIndices(mol);
                    // RDKit❗✔️:   auto volContribs = labuteContribs(mol);
                    // RDKit❗✔️:   std::vector<double> res(bins->size() + 1, 0.0);
                    // RDKit❗✔️:   for (unsigned int i = 0; i < propContribs.size(); ++i) {
                    // RDKit❗✔️:     auto nbin = std::upper_bound(bins->begin(), bins->end(), volContribs[i]) -
                    // RDKit❗✔️:                 bins->begin();
                    // RDKit❗✔️:     res[nbin] += propContribs[i];
                    // RDKit❗✔️:   }
                    // RDKit❗✔️:   return res;
                    // RDKit❗✔️: }
                    crate::assign_contribs_to_bins(&indices, &areas, AREA_CUTOFFS)?
                }
            };
            *self.slot(state) = Some(result);
        }
        Ok(self
            .slot(state)
            .as_deref()
            .expect("completed EState projection"))
    }

    fn slot(self, state: &mut crate::DescriptorComputedState) -> &mut Option<Vec<f64>> {
        match self {
            Self::AreaByEState => &mut state.area_by_estate,
            Self::EStateByArea => &mut state.estate_by_area,
        }
    }
}
