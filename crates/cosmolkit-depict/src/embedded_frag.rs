//! Private RDKit-style incremental 2D fragment state.

use std::cmp::Ordering;
use std::collections::BTreeMap;
use std::f64::consts::PI;
use std::sync::Mutex;

use cosmolkit_core::{
    DenseMatrix, MatrixError, PathError, RingInfo, TopologicalDistanceMatrixParams, shortest_path,
    topological_distance_matrix,
};
use cosmolkit_model::{AtomId, BondId, BondOrder, BondStereo, Hybridization, TopologyBlock};
use cosmolkit_search::{
    SearchTarget, SubstructMatchError, SubstructMatchParams, try_get_substruct_matches_with_params,
};

use crate::templates::{CoordinateTemplates, RingTemplate};

use crate::geometry::{
    BOND_LEN, GeometryError, Point2, PointMap, Transform2D, compute_bisect_point, embed_ring,
    rank_atoms_by_rank, reflect_point, set_neighbor_order,
};

// RDGeneral/utils.cpp owns one process-global generator. The mutex preserves
// that state without introducing a Rust data race; a positive seed resets it
// at the source call site, while zero/negative seeds continue its stream.
pub(crate) static DEPICT_RANDOM: Mutex<DepictRandom> = Mutex::new(DepictRandom { state: 42 });

#[derive(Debug, Clone, Copy)]
pub(crate) struct DepictRandom {
    state: u32,
}

impl DepictRandom {
    pub(crate) fn with_seed(seed: u32) -> Self {
        let mut generator = Self { state: 1 };
        generator.seed(seed);
        generator
    }

    pub(crate) fn seed(&mut self, seed: u32) {
        // Boost❗✔️: if(modulus == 0) {
        // Boost❗✔️:     _x = x0;
        // Boost❗✔️: } else {
        // Boost❗✔️:     _x = x0 % modulus;
        // Boost❗✔️: }
        // Boost❗✔️: if(increment == 0 && _x == 0) {
        // Boost❗✔️:     _x = 1;
        // Boost❗✔️: }
        // Behavior: minstd_rand has modulus 2147483647 and increment zero.
        // Complexity: constant-time integer remainder, as in Boost.
        self.state = seed % 2_147_483_647;
        if self.state == 0 {
            self.state = 1;
        }
    }

    pub(crate) fn reset_if_positive(&mut self, seed: i32) {
        // RDKit❗✔️: rng_type &getRandomGenerator(int seed) {
        // RDKit❗✔️:   if (seed > 0) {
        // RDKit❗✔️:     generator.seed(seed);
        // RDKit❗✔️:   }
        // RDKit❗✔️:   return generator;
        // RDKit❗✔️: }
        // Behavior: zero and negative seeds leave the shared stream untouched.
        // Complexity: constant-time check and at most one seed reduction.
        if seed > 0 {
            self.seed(seed as u32);
        }
    }

    pub(crate) fn next_raw(&mut self) -> u32 {
        // Boost❗✔️: typedef linear_congruential_engine<uint32_t, 48271, 0, 2147483647> minstd_rand;
        // Boost❗✔️: IntType operator()()
        // Boost❗✔️: {
        // Boost❗✔️:     _x = const_mod<IntType, m>::mult_add(a, _x, c);
        // Boost❗✔️:     return _x;
        // Boost❗✔️: }
        // Behavior: u64 intermediate exactly covers 48271 * (modulus - 1).
        // Complexity: one multiplication and remainder per source draw.
        self.state = ((u64::from(self.state) * 48_271) % 2_147_483_647) as u32;
        self.state
    }

    pub(crate) fn next_index(&mut self, count: usize) -> usize {
        assert!(count > 0 && count <= i32::MAX as usize);
        self.uniform_inclusive((count - 1) as u64) as usize
    }

    fn uniform_inclusive(&mut self, range: u64) -> u64 {
        // Boost❗✔️: const range_type range = random::detail::subtract<result_type>()(max_value, min_value);
        // Boost❗✔️: const base_result bmin = (eng.min)();
        // Boost❗✔️: const base_unsigned brange =
        // Boost❗✔️:   random::detail::subtract<base_result>()((eng.max)(), (eng.min)());
        // Boost❗✔️: if(range == 0) {
        // Boost❗✔️:   return min_value;
        // Boost❗✔️: } else if(brange == range) {
        // Boost❗✔️:   base_unsigned v = random::detail::subtract<base_result>()(eng(), bmin);
        // Boost❗✔️:   return random::detail::add<base_unsigned, result_type>()(v, min_value);
        // Boost❗✔️: } else if(brange < range) {
        // Boost❗✔️:   for(;;) {
        // Boost❗✔️:     range_type limit = (range+1)/(range_type(brange)+1);
        // Boost❗✔️:     range_type result = range_type(0);
        // Boost❗✔️:     range_type mult = range_type(1);
        // Boost❗✔️:     while(mult <= limit) {
        // Boost❗✔️:       result += static_cast<range_type>(static_cast<range_type>(random::detail::subtract<base_result>()(eng(), bmin)) * mult);
        // Boost❗✔️:       if(mult * range_type(brange) == range - mult + 1) return(result);
        // Boost❗✔️:       mult *= range_type(brange)+range_type(1);
        // Boost❗✔️:     }
        // Boost❗✔️:     range_type result_increment = generate_uniform_int(eng, 0, range/mult, boost::true_type());
        // Boost❗✔️:     result_increment *= mult;
        // Boost❗✔️:     result += result_increment;
        // Boost❗✔️:     if(result > range) continue;
        // Boost❗✔️:     return result;
        // Boost❗✔️:   }
        // Boost❗✔️: } else {
        // Boost❗✔️:   bucket_size = static_cast<mixed_range_type>(brange + 1) / (static_cast<mixed_range_type>(range)+1);
        // Boost❗✔️:   for(;;) {
        // Boost❗✔️:     mixed_range_type result = random::detail::subtract<base_result>()(eng(), bmin);
        // Boost❗✔️:     result /= bucket_size;
        // Boost❗✔️:     if(result <= static_cast<mixed_range_type>(range)) return result;
        // Boost❗✔️:   }
        // Boost❗✔️: }
        // Behavior: all reachable signed-int destination ranges are at most
        // i32::MAX-1, so u64 intermediates cannot overflow. The source's
        // bounded-overflow guards in its wider generic branch are unreachable.
        // Complexity: source rejection loops and draw count, no modulo bias.
        const BRANGE: u64 = 2_147_483_645;
        if range == 0 {
            return 0;
        }
        if range == BRANGE {
            return u64::from(self.next_raw() - 1);
        }
        if range > BRANGE {
            loop {
                let limit = (range + 1) / (BRANGE + 1);
                let mut result = 0;
                let mut mult = 1;
                while mult <= limit {
                    result += u64::from(self.next_raw() - 1) * mult;
                    if mult * BRANGE == range - mult + 1 {
                        return result;
                    }
                    mult *= BRANGE + 1;
                }
                result += self.uniform_inclusive(range / mult) * mult;
                if result <= range {
                    return result;
                }
            }
        }
        let bucket_size = (BRANGE + 1) / (range + 1);
        loop {
            let result = u64::from(self.next_raw() - 1) / bucket_size;
            if result <= range {
                return result;
            }
        }
    }

    pub(crate) fn next_real(&mut self) -> f64 {
        // Boost❗✔️: for(;;) {
        // Boost❗✔️:     result_type numerator = static_cast<T>(subtract<base_result>()(eng(), (eng.min)()));
        // Boost❗✔️:     result_type divisor = static_cast<T>(subtract<base_result>()((eng.max)(), (eng.min)())) + 1;
        // Boost❗✔️:     T result = numerator / divisor * (max_value - min_value) + min_value;
        // Boost❗✔️:     if(result < max_value) return result;
        // Boost❗✔️: }
        // Behavior: the fixed distribution is [0,1); draw again if rounding
        // reaches the upper endpoint. Complexity: one draw normally.
        loop {
            let value = f64::from(self.next_raw() - 1) / 2_147_483_646.0;
            if value < 1.0 {
                return value;
            }
        }
    }
}

#[derive(Debug, Clone, PartialEq)]
pub(crate) enum FragmentError {
    AtomIndexOutOfRange { atom: usize, atom_count: usize },
    AtomAlreadyEmbedded { atom: usize },
    AtomNotEmbedded { atom: usize },
    NotEnoughEmbeddedNeighbors { atom: usize, count: usize },
    CoincidentPoints,
    InvalidAngle,
    NoCommonAtoms,
    MismatchedTopology,
    EmptyAttachment { atom: usize },
    NonTetrahedralNoLigand { centre: usize },
    NonTetrahedralLigandOverflow { centre: usize },
    CisTransBondInvalid { bond: usize },
    CollisionBondInvalid { bond: usize },
    UndefinedSamplingDistance { first: usize, second: usize },
    Geometry(GeometryError),
    TemplateMatch(SubstructMatchError),
    GraphPath(PathError),
    GraphDistance(MatrixError),
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct DegreeFourCandidate {
    centre: usize,
    neighbor_atoms: Vec<usize>,
    neighbor_bonds: Vec<usize>,
}

impl From<GeometryError> for FragmentError {
    fn from(value: GeometryError) -> Self {
        Self::Geometry(value)
    }
}

impl From<SubstructMatchError> for FragmentError {
    fn from(value: SubstructMatchError) -> Self {
        Self::TemplateMatch(value)
    }
}

impl From<PathError> for FragmentError {
    fn from(value: PathError) -> Self {
        Self::GraphPath(value)
    }
}

impl From<MatrixError> for FragmentError {
    fn from(value: MatrixError) -> Self {
        Self::GraphDistance(value)
    }
}

fn ring_system_degree_counts(topology: &TopologyBlock, included: &[bool]) -> [usize; 5] {
    // BEGIN RDKIT CPP FUNCTION EmbeddedFrag::matchToTemplate degreeCounts
    // RDKit❗❗:       std::array<int, 5> degrees_count({0, 0, 0, 0, 0});
    // RDKit❗❗:       for (auto atom : mol.atoms()) {
    // RDKit❗❗:         if (atom->getAtomicNum() == DUMMY_ATOMIC_NUM) {
    // RDKit❗❗:           continue;
    // RDKit❗❗:         }
    // RDKit❗❗:         auto degree = 0u;
    // RDKit❗❗:         for (auto nbr : mol.atomNeighbors(atom)) {
    // RDKit❗❗:           if (nbr->getAtomicNum() != DUMMY_ATOMIC_NUM) {
    // RDKit❗❗:             ++degree;
    // RDKit❗❗:             if (degree == 4) {
    // RDKit❗❗:               break;
    // RDKit❗❗:             }
    // RDKit❗❗:           }
    // RDKit❗❗:         }
    // RDKit❗❗:         degrees_count[degree]++;
    // RDKit❗❗:       }
    // RDKit❗❗:       return degrees_count;
    // RDKit❗❗:     };
    // RDKit❗❗:     if (degreeCounts(rs_mol) != degreeCounts(*mol)) {
    // END RDKIT CPP FUNCTION EmbeddedFrag::matchToTemplate degreeCounts
    // Behavior: counts only non-sentinel neighbors and caps at degree four.
    // Complexity: one adjacency scan per participating atom, as in source.
    let mut counts = [0; 5];
    for (index, &inside) in included.iter().enumerate() {
        if !inside {
            continue;
        }
        let degree = topology
            .adjacency
            .neighbors_of(index)
            .iter()
            .filter(|neighbor| included[neighbor.atom_index])
            .take(4)
            .count();
        counts[degree] += 1;
    }
    counts
}

fn template_degree_counts(template: &RingTemplate) -> [usize; 5] {
    // BEGIN RDKIT CPP FUNCTION EmbeddedFrag::matchToTemplate degreeCounts
    // RDKit❗❗:       std::array<int, 5> degrees_count({0, 0, 0, 0, 0});
    // RDKit❗❗:       for (auto atom : mol.atoms()) {
    // RDKit❗❗:         if (atom->getAtomicNum() == DUMMY_ATOMIC_NUM) {
    // RDKit❗❗:           continue;
    // RDKit❗❗:         }
    // RDKit❗❗:         auto degree = 0u;
    // RDKit❗❗:         for (auto nbr : mol.atomNeighbors(atom)) {
    // RDKit❗❗:           if (nbr->getAtomicNum() != DUMMY_ATOMIC_NUM) {
    // RDKit❗❗:             ++degree;
    // RDKit❗❗:             if (degree == 4) {
    // RDKit❗❗:               break;
    // RDKit❗❗:             }
    // RDKit❗❗:           }
    // RDKit❗❗:         }
    // RDKit❗❗:         degrees_count[degree]++;
    // RDKit❗❗:       }
    // RDKit❗❗:       return degrees_count;
    // RDKit❗❗:     };
    // RDKit❗❗:     if (degreeCounts(rs_mol) != degreeCounts(*mol)) {
    // END RDKIT CPP FUNCTION EmbeddedFrag::matchToTemplate degreeCounts
    // Behavior: counts only non-sentinel neighbors and caps at degree four.
    // Complexity: one adjacency scan per participating atom, as in source.
    let mut counts = [0; 5];
    for (index, atom) in template.query.atoms().iter().enumerate() {
        if atom.atomic_number() == 200 {
            continue;
        }
        let degree = template.query.adjacency()[index]
            .iter()
            .filter(|&&(neighbor, _)| template.query.atoms()[neighbor].atomic_number() != 200)
            .take(4)
            .count();
        counts[degree] += 1;
    }
    counts
}

fn template_point(template: &RingTemplate, index: usize) -> Point2 {
    // RDKit❗❗: const auto &conf = template_mol->getConformer();
    // RDKit❗❗: for (auto &[template_aidx, rs_aidx] : match) {
    // RDKit❗❗:   EmbeddedAtom new_at(rs_aidx, conf.getAtomPos(template_aidx));
    // RDKit❗❗:   new_at.df_fixed = true;
    // RDKit❗❗:   d_eatoms.emplace(rs_aidx, new_at);
    // RDKit❗❗: }
    // Behavior: projection selects the first stored template conformer.
    // Complexity: O(1) indexed coordinate access without allocation.
    if let Some(rows) = template.query.coordinates_2d() {
        rows[index]
    } else {
        let row = template.query.conformers_3d()[0].coordinates()[index];
        [row[0], row[1]]
    }
}

fn check_template_stereo(
    topology: &TopologyBlock,
    template: &RingTemplate,
    mapping: &[usize],
) -> bool {
    // BEGIN RDKIT CPP FUNCTION checkStereoChemistry
    // RDKit❗❗: static bool checkStereoChemistry(const RDKit::ROMol &mol,
    // RDKit❗❗:                                  const RDKit::ROMol &template_mol,
    // RDKit❗❗:                                  RDKit::MatchVectType match) {
    // RDKit❗❗:   for (auto bond : mol.bonds()) {
    // RDKit❗❗:     if (bond->getBondType() != RDKit::Bond::DOUBLE ||
    // RDKit❗❗:         bond->getStereo() == RDKit::Bond::STEREOANY ||
    // RDKit❗❗:         bond->getStereo() == RDKit::Bond::STEREONONE) {
    // RDKit❗❗:       continue;
    // RDKit❗❗:     }
    // RDKit❗❗:     // get the four atoms around the double bond
    // RDKit❗❗:     auto neighbors = bond->getStereoAtoms();
    // RDKit❗❗:     if (neighbors.size() != 2) {
    // RDKit❗❗:       continue;
    // RDKit❗❗:     }
    // RDKit❗❗:     int atom1_neighbor1 = neighbors[0];
    // RDKit❗❗:     int atom2_neighbor1 = neighbors[1];
    // RDKit❗❗:     int atom1 = bond->getBeginAtomIdx();
    // RDKit❗❗:     int atom2 = bond->getEndAtomIdx();
    // RDKit❗❗:
    // RDKit❗❗:     // now get the other two atoms that are not part of the double bond (if any)
    // RDKit❗❗:     int atom1_neighbor2 = -1;
    // RDKit❗❗:     int atom2_neighbor2 = -1;
    // RDKit❗❗:     if (mol.getAtomWithIdx(atom1)->getDegree() > 2) {
    // RDKit❗❗:       for (auto neighbor : mol.atomNeighbors(mol.getAtomWithIdx(atom1))) {
    // RDKit❗❗:         if (static_cast<int>(neighbor->getIdx()) != atom1_neighbor1 &&
    // RDKit❗❗:             static_cast<int>(neighbor->getIdx()) != atom2) {
    // RDKit❗❗:           atom1_neighbor2 = neighbor->getIdx();
    // RDKit❗❗:           break;
    // RDKit❗❗:         }
    // RDKit❗❗:       }
    // RDKit❗❗:     }
    // RDKit❗❗:     if (mol.getAtomWithIdx(atom2)->getDegree() > 2) {
    // RDKit❗❗:       for (auto neighbor : mol.atomNeighbors(mol.getAtomWithIdx(atom2))) {
    // RDKit❗❗:         if (static_cast<int>(neighbor->getIdx()) != atom2_neighbor1 &&
    // RDKit❗❗:             static_cast<int>(neighbor->getIdx()) != atom1) {
    // RDKit❗❗:           atom2_neighbor2 = neighbor->getIdx();
    // RDKit❗❗:           break;
    // RDKit❗❗:         }
    // RDKit❗❗:       }
    // RDKit❗❗:     }
    // RDKit❗❗:
    // RDKit❗❗:     // find the template atoms that correspond to the four atoms
    // RDKit❗❗:     int template_atom1 = -1;
    // RDKit❗❗:     int template_atom2 = -1;
    // RDKit❗❗:     int template_atom1_neighbor1 = -1;
    // RDKit❗❗:     int template_atom1_neighbor2 = -1;
    // RDKit❗❗:     int template_atom2_neighbor1 = -1;
    // RDKit❗❗:     int template_atom2_neighbor2 = -1;
    // RDKit❗❗:     for (auto &[template_aidx, rs_aidx] : match) {
    // RDKit❗❗:       if (rs_aidx == atom1) {
    // RDKit❗❗:         template_atom1 = template_aidx;
    // RDKit❗❗:       } else if (rs_aidx == atom2) {
    // RDKit❗❗:         template_atom2 = template_aidx;
    // RDKit❗❗:       } else if (rs_aidx == atom1_neighbor1) {
    // RDKit❗❗:         template_atom1_neighbor1 = template_aidx;
    // RDKit❗❗:       } else if (rs_aidx == atom2_neighbor1) {
    // RDKit❗❗:         template_atom2_neighbor1 = template_aidx;
    // RDKit❗❗:       } else if (rs_aidx == atom1_neighbor2) {
    // RDKit❗❗:         template_atom1_neighbor2 = template_aidx;
    // RDKit❗❗:       } else if (rs_aidx == atom2_neighbor2) {
    // RDKit❗❗:         template_atom2_neighbor2 = template_aidx;
    // RDKit❗❗:       }
    // RDKit❗❗:     }
    // RDKit❗❗:
    // RDKit❗❗:     // there's a chance that the atoms controlling the double bond stereochem in
    // RDKit❗❗:     // the molecule are not the atoms that matched to the template, handle that
    // RDKit❗❗:     // here by swapping to the other atom
    // RDKit❗❗:     bool swapStereo = false;
    // RDKit❗❗:     if (template_atom1_neighbor1 == -1) {
    // RDKit❗❗:       template_atom1_neighbor1 = template_atom1_neighbor2;
    // RDKit❗❗:       swapStereo = !swapStereo;
    // RDKit❗❗:     }
    // RDKit❗❗:     if (template_atom2_neighbor1 == -1) {
    // RDKit❗❗:       template_atom2_neighbor1 = template_atom2_neighbor2;
    // RDKit❗❗:       swapStereo = !swapStereo;
    // RDKit❗❗:     }
    // RDKit❗❗:
    // RDKit❗❗:     if (template_atom1 == -1 || template_atom2 == -1 ||
    // RDKit❗❗:         template_atom1_neighbor1 == -1 || template_atom2_neighbor1 == -1) {
    // RDKit❗❗:       return false;
    // RDKit❗❗:     }
    // RDKit❗❗:
    // RDKit❗❗:     const auto &conf = template_mol.getConformer();
    // RDKit❗❗:     const auto &atom1_loc = conf.getAtomPos(template_atom1);
    // RDKit❗❗:     const auto &atom2_loc = conf.getAtomPos(template_atom2);
    // RDKit❗❗:     const auto &atom1_neighbor_loc = conf.getAtomPos(template_atom1_neighbor1);
    // RDKit❗❗:     const auto &atom2_neighbor_loc = conf.getAtomPos(template_atom2_neighbor1);
    // RDKit❗❗:     // check if the two neighbors are on the same side of the bond
    // RDKit❗❗:     const auto v12 = atom1_neighbor_loc - atom1_loc;
    // RDKit❗❗:     const auto v42 = atom2_neighbor_loc - atom1_loc;
    // RDKit❗❗:     const auto v32 = atom2_loc - atom1_loc;
    // RDKit❗❗:     auto cross1 = v32.x * v12.y - v32.y * v12.x;
    // RDKit❗❗:     auto cross2 = v32.x * v42.y - v32.y * v42.x;
    // RDKit❗❗:     bool is_cis = cross1 * cross2 > 0;
    // RDKit❗❗:     if (swapStereo) {
    // RDKit❗❗:       is_cis = !is_cis;
    // RDKit❗❗:     }
    // RDKit❗❗:     if (is_cis != (bond->getStereo() == RDKit::Bond::STEREOZ ||
    // RDKit❗❗:                    bond->getStereo() == RDKit::Bond::STEREOCIS)) {
    // RDKit❗❗:       return false;
    // RDKit❗❗:     }
    // RDKit❗❗:   }
    // RDKit❗❗:   return true;
    // RDKit❗❗: }
    // END RDKIT CPP FUNCTION checkStereoChemistry
    // Behavior: full source mapping/swap/cross-product order is ported; the
    // regression gate must establish parity for the supported stereo states.
    // Complexity: both implementations scan bonds and the match vector per
    // stereo bond; this function allocates no additional per-bond collection.
    for bond in &topology.bonds {
        if bond.order() != cosmolkit_model::BondOrder::Double
            || matches!(
                bond.stereo(),
                cosmolkit_model::BondStereo::Any | cosmolkit_model::BondStereo::None
            )
        {
            continue;
        }
        let Some([neighbor1, neighbor2]) = bond.stereo_atoms() else {
            continue;
        };
        let atom1 = bond.begin().index();
        let atom2 = bond.end().index();
        let other1 = (topology.adjacency.neighbors_of(atom1).len() > 2)
            .then(|| {
                topology
                    .adjacency
                    .neighbors_of(atom1)
                    .iter()
                    .find_map(|row| {
                        (row.atom_index != neighbor1.index() && row.atom_index != atom2)
                            .then_some(row.atom_index)
                    })
            })
            .flatten();
        let other2 = (topology.adjacency.neighbors_of(atom2).len() > 2)
            .then(|| {
                topology
                    .adjacency
                    .neighbors_of(atom2)
                    .iter()
                    .find_map(|row| {
                        (row.atom_index != neighbor2.index() && row.atom_index != atom1)
                            .then_some(row.atom_index)
                    })
            })
            .flatten();
        let find = |target: usize| mapping.iter().position(|&mapped| mapped == target);
        let Some(template_atom1) = find(atom1) else {
            return false;
        };
        let Some(template_atom2) = find(atom2) else {
            return false;
        };
        let mut template_neighbor1 = find(neighbor1.index());
        let mut template_neighbor2 = find(neighbor2.index());
        let mut swap = false;
        if template_neighbor1.is_none() {
            template_neighbor1 = other1.and_then(find);
            swap = !swap;
        }
        if template_neighbor2.is_none() {
            template_neighbor2 = other2.and_then(find);
            swap = !swap;
        }
        let (Some(template_neighbor1), Some(template_neighbor2)) =
            (template_neighbor1, template_neighbor2)
        else {
            return false;
        };
        let origin = template_point(template, template_atom1);
        let end = template_point(template, template_atom2);
        let reference1 = template_point(template, template_neighbor1);
        let reference2 = template_point(template, template_neighbor2);
        let axis = [end[0] - origin[0], end[1] - origin[1]];
        let arm1 = [reference1[0] - origin[0], reference1[1] - origin[1]];
        let arm2 = [reference2[0] - origin[0], reference2[1] - origin[1]];
        let cross1 = axis[0] * arm1[1] - axis[1] * arm1[0];
        let cross2 = axis[0] * arm2[1] - axis[1] * arm2[0];
        let mut cis = cross1 * cross2 > 0.0;
        if swap {
            cis = !cis;
        }
        if cis
            != matches!(
                bond.stereo(),
                cosmolkit_model::BondStereo::Z | cosmolkit_model::BondStereo::Cis
            )
        {
            return false;
        }
    }
    true
}

fn ring_union(rings: &[Vec<usize>]) -> Vec<usize> {
    // RDKit❗✔️: void Union(const VECT_INT_VECT &rings, INT_VECT &res, const INT_VECT *exclude) {
    // RDKit❗✔️:   res.resize(0);
    // RDKit❗✔️:   INT_VECT ring;
    // RDKit❗✔️:   unsigned int id;
    // RDKit❗✔️:   auto nrings = static_cast<unsigned int>(rings.size());
    // RDKit❗✔️:   INT_VECT_CI ri;
    // RDKit❗✔️:   for (id = 0; id < nrings; id++) {
    // RDKit❗✔️:     if (exclude) {
    // RDKit❗✔️:       if (std::find(exclude->begin(), exclude->end(), static_cast<int>(id)) !=
    // RDKit❗✔️:           exclude->end()) {
    // RDKit❗✔️:         continue;
    // RDKit❗✔️:       }
    // RDKit❗✔️:     }
    // RDKit❗✔️:     ring = rings[id];
    // RDKit❗✔️:     for (ri = ring.begin(); ri != ring.end(); ri++) {
    // RDKit❗✔️:       if (std::find(res.begin(), res.end(), (*ri)) == res.end()) {
    // RDKit❗✔️:         res.push_back(*ri);
    // RDKit❗✔️:       }
    // RDKit❗✔️:     }
    // RDKit❗✔️:   }
    // RDKit❗✔️: }
    // Behavior: source-order first occurrence, with no sorting.
    // Complexity: linear membership scans retain the source's quadratic bound.
    let mut result = Vec::new();
    for ring in rings {
        for &atom in ring {
            if !result.contains(&atom) {
                result.push(atom);
            }
        }
    }
    result
}

pub(crate) fn ring_systems(rings: &[Vec<usize>]) -> Vec<Vec<usize>> {
    // RDKit❗✔️: void makeRingNeighborMap(const VECT_INT_VECT &brings,
    // RDKit❗✔️:                          INT_INT_VECT_MAP &neighMap, unsigned int maxSize,
    // RDKit❗✔️:                          unsigned int maxOverlapSize) {
    // RDKit❗✔️:   auto nrings = rdcast<int>(brings.size());
    // RDKit❗✔️:   int i, j;
    // RDKit❗✔️:   INT_VECT ring1;
    // RDKit❗✔️:   for (i = 0; i < nrings; ++i) {
    // RDKit❗✔️:     neighMap[i];
    // RDKit❗✔️:     if (maxSize && brings[i].size() > maxSize) {
    // RDKit❗✔️:       continue;
    // RDKit❗✔️:     }
    // RDKit❗✔️:     ring1 = brings[i];
    // RDKit❗✔️:     for (j = i + 1; j < nrings; ++j) {
    // RDKit❗✔️:       if (maxSize && brings[j].size() > maxSize) {
    // RDKit❗✔️:         continue;
    // RDKit❗✔️:       }
    // RDKit❗✔️:       INT_VECT inter;
    // RDKit❗✔️:       Intersect(ring1, brings[j], inter);
    // RDKit❗✔️:       if (inter.size() > 0 &&
    // RDKit❗✔️:           (!maxOverlapSize || inter.size() <= maxOverlapSize)) {
    // RDKit❗✔️:         neighMap[i].push_back(j);
    // RDKit❗✔️:         neighMap[j].push_back(i);
    // RDKit❗✔️:       }
    // RDKit❗✔️:     }
    // RDKit❗✔️:   }
    // RDKit❗✔️: }
    // RDKit❗✔️: void pickFusedRings(int curr, const INT_INT_VECT_MAP &neighMap, INT_VECT &res,
    // RDKit❗✔️:                     boost::dynamic_bitset<> &done, int depth) {
    // RDKit❗✔️:   auto pos = neighMap.find(curr);
    // RDKit❗✔️:   PRECONDITION(pos != neighMap.end(), "bad argument");
    // RDKit❗✔️:   done[curr] = 1;
    // RDKit❗✔️:   res.push_back(curr);
    // RDKit❗✔️:   const auto &neighs = pos->second;
    // RDKit❗✔️:   for (int neigh : neighs) {
    // RDKit❗✔️:     if (!done[neigh]) {
    // RDKit❗✔️:       pickFusedRings(neigh, neighMap, res, done, depth + 1);
    // RDKit❗✔️:     }
    // RDKit❗✔️:   }
    // RDKit❗✔️: }
    // Behavior: no max-size/overlap filtering at this call site; spiro rings
    // share an atom and therefore belong to one embedding system.
    // Complexity: pairwise ring scans and one traversal match the source;
    // explicit DFS stack avoids unbounded Rust call-stack growth.
    let mut neighbors = vec![Vec::new(); rings.len()];
    for left in 0..rings.len() {
        for right in left + 1..rings.len() {
            if rings[left].iter().any(|atom| rings[right].contains(atom)) {
                neighbors[left].push(right);
                neighbors[right].push(left);
            }
        }
    }
    let mut done = vec![false; rings.len()];
    let mut systems = Vec::new();
    for root in 0..rings.len() {
        if done[root] {
            continue;
        }
        let mut system = Vec::new();
        let mut stack = vec![(root, 0)];
        done[root] = true;
        system.push(root);
        while let Some((ring, next)) = stack.last_mut() {
            if *next == neighbors[*ring].len() {
                stack.pop();
                continue;
            }
            let adjacent = neighbors[*ring][*next];
            *next += 1;
            if !done[adjacent] {
                done[adjacent] = true;
                system.push(adjacent);
                stack.push((adjacent, 0));
            }
        }
        systems.push(system);
    }
    systems
}

fn pick_first_ring_to_embed(topology: &TopologyBlock, rings: &[Vec<usize>]) -> usize {
    // RDKit❗✔️: int pickFirstRingToEmbed(const RDKit::ROMol &mol,
    // RDKit❗✔️:                          const RDKit::VECT_INT_VECT &fusedRings) {
    // RDKit❗✔️:   int res = -1;
    // RDKit❗✔️:   unsigned int maxSize = 0;
    // RDKit❗✔️:   int subs, minsubs = static_cast<int>(1e8);
    // RDKit❗✔️:   int cnt = 0;
    // RDKit❗✔️:   for (const auto &fusedRing : fusedRings) {
    // RDKit❗✔️:     subs = 0;
    // RDKit❗✔️:     for (auto rii : fusedRing) {
    // RDKit❗✔️:       if (mol.getAtomWithIdx(rii)->getDegree() > 2) {
    // RDKit❗✔️:         ++subs;
    // RDKit❗✔️:       }
    // RDKit❗✔️:     }
    // RDKit❗✔️:     if (subs < minsubs) {
    // RDKit❗✔️:       res = cnt;
    // RDKit❗✔️:       minsubs = subs;
    // RDKit❗✔️:       maxSize = fusedRing.size();
    // RDKit❗✔️:     } else if (subs == minsubs) {
    // RDKit❗✔️:       if (fusedRing.size() > maxSize) {
    // RDKit❗✔️:         res = cnt;
    // RDKit❗✔️:         maxSize = fusedRing.size();
    // RDKit❗✔️:       }
    // RDKit❗✔️:     }
    // RDKit❗✔️:     cnt++;
    // RDKit❗✔️:   }
    // RDKit❗✔️:   return res;
    // RDKit❗✔️: }
    // Behavior: minimum substituted atoms, then largest ring, then first.
    // Complexity: one adjacency degree lookup per ring atom.
    let mut best = 0;
    let mut min_substituents = usize::MAX;
    let mut max_size = 0;
    for (index, ring) in rings.iter().enumerate() {
        let substituents = ring
            .iter()
            .filter(|&&atom| topology.adjacency.neighbors_of(atom).len() > 2)
            .count();
        if substituents < min_substituents
            || (substituents == min_substituents && ring.len() > max_size)
        {
            best = index;
            min_substituents = substituents;
            max_size = ring.len();
        }
    }
    best
}

fn find_core_rings(topology: &TopologyBlock, rings: &[Vec<usize>]) -> Vec<usize> {
    // RDKit❗✔️: RDKit::VECT_INT_VECT findCoreRings(const RDKit::VECT_INT_VECT &fusedRings,
    // RDKit❗✔️:                                    RDKit::INT_VECT &coreRingsIds,
    // RDKit❗✔️:                                    const RDKit::ROMol &mol) {
    // RDKit❗✔️:   boost::dynamic_bitset<> removedRings(fusedRings.size());
    // RDKit❗✔️:   bool removedARing = false;
    // RDKit❗✔️:   do {
    // RDKit❗✔️:     removedARing = false;
    // RDKit❗✔️:     for (unsigned int currRingId = 0; currRingId < fusedRings.size();
    // RDKit❗✔️:          currRingId++) {
    // RDKit❗✔️:       if (removedRings[currRingId] || removedARing) {
    // RDKit❗✔️:         continue;
    // RDKit❗✔️:       }
    // RDKit❗✔️:       auto nIntersectingAtoms = 0u;
    // RDKit❗✔️:       int aid1 = -1;
    // RDKit❗✔️:       int aid2 = -1;
    // RDKit❗✔️:       for (unsigned int otherRingId = 0; otherRingId < fusedRings.size();
    // RDKit❗✔️:            otherRingId++) {
    // RDKit❗✔️:         if (currRingId == otherRingId || removedRings[otherRingId]) {
    // RDKit❗✔️:           continue;
    // RDKit❗✔️:         }
    // RDKit❗✔️:         RDKit::INT_VECT commmonAtoms;
    // RDKit❗✔️:         RDKit::Intersect(fusedRings[currRingId], fusedRings[otherRingId],
    // RDKit❗✔️:                          commmonAtoms);
    // RDKit❗✔️:         for (auto rii : commmonAtoms) {
    // RDKit❗✔️:           if (rii != aid1 && rii != aid2) {
    // RDKit❗✔️:             ++nIntersectingAtoms;
    // RDKit❗✔️:             if (aid1 == -1) {
    // RDKit❗✔️:               aid1 = rii;
    // RDKit❗✔️:             } else {
    // RDKit❗✔️:               aid2 = rii;
    // RDKit❗✔️:             }
    // RDKit❗✔️:             if (nIntersectingAtoms == 2) {
    // RDKit❗✔️:               break;
    // RDKit❗✔️:             }
    // RDKit❗✔️:           }
    // RDKit❗✔️:         }
    // RDKit❗✔️:       }
    // RDKit❗✔️:       if (nIntersectingAtoms == 1 ||
    // RDKit❗✔️:           (nIntersectingAtoms == 2 &&
    // RDKit❗✔️:            mol.getBondBetweenAtoms(aid1, aid2) != nullptr)) {
    // RDKit❗✔️:         removedRings[currRingId] = true;
    // RDKit❗✔️:         removedARing = true;
    // RDKit❗✔️:       }
    // RDKit❗✔️:     }
    // RDKit❗✔️:   } while (removedARing);
    // RDKit❗✔️:   RDKit::VECT_INT_VECT res;
    // RDKit❗✔️:   for (unsigned int currRingId = 0; currRingId < fusedRings.size();
    // RDKit❗✔️:        currRingId++) {
    // RDKit❗✔️:     if (!removedRings[currRingId]) {
    // RDKit❗✔️:       res.push_back(fusedRings[currRingId]);
    // RDKit❗✔️:       coreRingsIds.push_back(currRingId);
    // RDKit❗✔️:     }
    // RDKit❗✔️:   }
    // RDKit❗✔️:   return res;
    // RDKit❗✔️: }
    // Behavior: one removal per iteration, including the source's unusual
    // intersection counter that can exceed two across multiple other rings.
    // Complexity: source-order repeated ring intersections; the bond lookup
    // uses adjacency rather than a full bond scan.
    let mut removed = vec![false; rings.len()];
    loop {
        let mut removed_one = false;
        for current in 0..rings.len() {
            if removed[current] || removed_one {
                continue;
            }
            let mut intersection_count = 0;
            let mut first = None;
            let mut second = None;
            for other in 0..rings.len() {
                if current == other || removed[other] {
                    continue;
                }
                for &atom in &rings[current] {
                    if rings[other].contains(&atom) && Some(atom) != first && Some(atom) != second {
                        intersection_count += 1;
                        if first.is_none() {
                            first = Some(atom);
                        } else {
                            second = Some(atom);
                        }
                        if intersection_count == 2 {
                            break;
                        }
                    }
                }
            }
            let adjacent = match (first, second) {
                (Some(a), Some(b)) => topology
                    .adjacency
                    .neighbors_of(a)
                    .iter()
                    .any(|row| row.atom_index == b),
                _ => false,
            };
            if intersection_count == 1 || (intersection_count == 2 && adjacent) {
                removed[current] = true;
                removed_one = true;
            }
        }
        if !removed_one {
            break;
        }
    }
    removed
        .iter()
        .enumerate()
        .filter_map(|(index, &is_removed)| (!is_removed).then_some(index))
        .collect()
}

fn find_next_ring_to_embed(done: &[usize], rings: &[Vec<usize>]) -> (usize, Vec<usize>) {
    // RDKit❗✔️: RDKit::INT_VECT findNextRingToEmbed(const RDKit::INT_VECT &doneRings,
    // RDKit❗✔️:                                     const RDKit::VECT_INT_VECT &fusedRings,
    // RDKit❗✔️:                                     int &nextId) {
    // RDKit❗✔️:   PRECONDITION(doneRings.size() > 0, "");
    // RDKit❗✔️:   PRECONDITION(fusedRings.size() > 1, "");
    // RDKit❗✔️:   RDKit::INT_VECT commonAtoms, res, doneAtoms, notDone;
    // RDKit❗✔️:   for (int i = 0; i < rdcast<int>(fusedRings.size()); i++) {
    // RDKit❗✔️:     if (std::find(doneRings.begin(), doneRings.end(), i) == doneRings.end()) {
    // RDKit❗✔️:       notDone.push_back(i);
    // RDKit❗✔️:     }
    // RDKit❗✔️:   }
    // RDKit❗✔️:   RDKit::Union(fusedRings, doneAtoms, &notDone);
    // RDKit❗✔️:   int maxCommonAtoms = 0;
    // RDKit❗✔️:   int currRingId = 0;
    // RDKit❗✔️:   for (const auto &fusedRing : fusedRings) {
    // RDKit❗✔️:     if (std::find(doneRings.begin(), doneRings.end(), currRingId) !=
    // RDKit❗✔️:         doneRings.end()) {
    // RDKit❗✔️:       currRingId++;
    // RDKit❗✔️:       continue;
    // RDKit❗✔️:     }
    // RDKit❗✔️:     commonAtoms.clear();
    // RDKit❗✔️:     int numCommonAtoms = 0;
    // RDKit❗✔️:     for (auto rii : fusedRing) {
    // RDKit❗✔️:       if (std::find(doneAtoms.begin(), doneAtoms.end(), (rii)) !=
    // RDKit❗✔️:           doneAtoms.end()) {
    // RDKit❗✔️:         commonAtoms.push_back(rii);
    // RDKit❗✔️:         numCommonAtoms++;
    // RDKit❗✔️:       }
    // RDKit❗✔️:     }
    // RDKit❗✔️:     if (numCommonAtoms == 2) {
    // RDKit❗✔️:       nextId = currRingId;
    // RDKit❗✔️:       return commonAtoms;
    // RDKit❗✔️:     }
    // RDKit❗✔️:     if (numCommonAtoms > maxCommonAtoms) {
    // RDKit❗✔️:       maxCommonAtoms = numCommonAtoms;
    // RDKit❗✔️:       nextId = currRingId;
    // RDKit❗✔️:       res = commonAtoms;
    // RDKit❗✔️:     }
    // RDKit❗✔️:     ++currRingId;
    // RDKit❗✔️:   }
    // RDKit❗✔️:   unsigned int cmnLst = 0;
    // RDKit❗✔️:   unsigned int nCmn = res.size();
    // RDKit❗✔️:   for (unsigned int i = 0; i < nCmn; i++) {
    // RDKit❗✔️:     if (res[i] == fusedRings[nextId][i]) {
    // RDKit❗✔️:       cmnLst++;
    // RDKit❗✔️:     } else {
    // RDKit❗✔️:       break;
    // RDKit❗✔️:     }
    // RDKit❗✔️:   }
    // RDKit❗✔️:   if ((cmnLst > 0) && (cmnLst < res.size())) {
    // RDKit❗✔️:     RDKit::INT_VECT tempV = res;
    // RDKit❗✔️:     for (unsigned int i = cmnLst; i < nCmn; i++) {
    // RDKit❗✔️:       res[i - cmnLst] = tempV[i];
    // RDKit❗✔️:     }
    // RDKit❗✔️:     unsigned int nMov = nCmn - cmnLst;
    // RDKit❗✔️:     for (unsigned int i = 0; i < cmnLst; i++) {
    // RDKit❗✔️:       res[nMov + i] = tempV[i];
    // RDKit❗✔️:     }
    // RDKit❗✔️:   }
    // RDKit❗✔️:   POSTCONDITION(res.size() > 0, "");
    // RDKit❗✔️:   return res;
    // RDKit❗✔️: }
    // Behavior: exact two-atom overlap wins immediately; otherwise first
    // largest overlap, with source's wrapped-chain rotation.
    // Complexity: source-order linear membership tests; no duplicate graph.
    let done_rings: Vec<_> = rings
        .iter()
        .enumerate()
        .filter_map(|(index, ring)| done.contains(&index).then_some(ring.clone()))
        .collect();
    let done_atoms = ring_union(&done_rings);
    let mut selected = None;
    let mut common = Vec::new();
    for (index, ring) in rings.iter().enumerate() {
        if done.contains(&index) {
            continue;
        }
        let overlap: Vec<_> = ring
            .iter()
            .copied()
            .filter(|atom| done_atoms.contains(atom))
            .collect();
        if overlap.len() == 2 {
            return (index, overlap);
        }
        if overlap.len() > common.len() {
            selected = Some(index);
            common = overlap;
        }
    }
    let next = selected.expect("connected ring system has an overlapping next ring");
    let prefix = common
        .iter()
        .zip(&rings[next])
        .take_while(|(a, b)| a == b)
        .count();
    if prefix > 0 && prefix < common.len() {
        common.rotate_left(prefix);
    }
    (next, common)
}

fn mirror_trans_ring_atoms(topology: &TopologyBlock, ring: &[usize], coords: &mut PointMap) {
    // RDKit❗✔️: static void mirrorTransRingAtoms(const RDKit::ROMol &mol,
    // RDKit❗✔️:                                  const RDKit::INT_VECT &ring,
    // RDKit❗✔️:                                  RDGeom::INT_POINT2D_MAP &coords) {
    // RDKit❗✔️:   RDKit::INT_VECT transRingAtoms;
    // RDKit❗✔️:   for (size_t i = 0; i < ring.size(); ++i) {
    // RDKit❗✔️:     const auto atom1 = ring[i];
    // RDKit❗✔️:     const auto atom2 = ring[(i + 1) % ring.size()];
    // RDKit❗✔️:     const auto bond = mol.getBondBetweenAtoms(atom1, atom2);
    // RDKit❗✔️:     if (bond->getBondType() != RDKit::Bond::DOUBLE) {
    // RDKit❗✔️:       continue;
    // RDKit❗✔️:     }
    // RDKit❗✔️:     const auto stype = bond->getStereo();
    // RDKit❗✔️:     if (stype <= RDKit::Bond::STEREOANY) {
    // RDKit❗✔️:       continue;
    // RDKit❗✔️:     }
    // RDKit❗✔️:     const auto &neighbors = bond->getStereoAtoms();
    // RDKit❗✔️:     if (neighbors.size() != 2) {
    // RDKit❗✔️:       continue;
    // RDKit❗✔️:     }
    // RDKit❗✔️:     const auto leftIsIn =
    // RDKit❗✔️:         std::find(ring.begin(), ring.end(), neighbors[0]) != ring.end();
    // RDKit❗✔️:     const auto rightIsIn =
    // RDKit❗✔️:         std::find(ring.begin(), ring.end(), neighbors[1]) != ring.end();
    // RDKit❗✔️:     bool isTrans = false;
    // RDKit❗✔️:     if (stype == RDKit::Bond::STEREOTRANS || stype == RDKit::Bond::STEREOE) {
    // RDKit❗✔️:       if (leftIsIn == rightIsIn) {
    // RDKit❗✔️:         isTrans = true;
    // RDKit❗✔️:       }
    // RDKit❗✔️:     } else if (leftIsIn != rightIsIn) {
    // RDKit❗✔️:       isTrans = true;
    // RDKit❗✔️:     }
    // RDKit❗✔️:     if (!isTrans) {
    // RDKit❗✔️:       continue;
    // RDKit❗✔️:     }
    // RDKit❗✔️:     const auto left = ring[(i + ring.size() - 1) % ring.size()];
    // RDKit❗✔️:     const auto right = atom2;
    // RDKit❗✔️:     const auto last = coords[left];
    // RDKit❗✔️:     const auto ref = coords[right];
    // RDKit❗✔️:     const auto interest = coords[atom1];
    // RDKit❗✔️:     const auto d = last - ref;
    // RDKit❗✔️:     const double a = (d.x * d.x - d.y * d.y) / d.dotProduct(d);
    // RDKit❗✔️:     const double b = 2 * d.x * d.y / d.dotProduct(d);
    // RDKit❗✔️:     const double x =
    // RDKit❗✔️:         a * (interest.x - ref.x) + b * (interest.y - ref.y) + ref.x;
    // RDKit❗✔️:     const double y =
    // RDKit❗✔️:         b * (interest.x - ref.x) - a * (interest.y - ref.y) + ref.y;
    // RDKit❗✔️:     coords[atom1] = RDGeom::Point2D(x, y);
    // RDKit❗✔️:   }
    // RDKit❗✔️: }
    // Behavior: only source-defined trans-with-respect-to-ring double bonds
    // are mirrored; the other cyclic coordinates stay bitwise untouched.
    // Complexity: adjacency lookup and one ring membership scan per edge.
    for (index, &atom1) in ring.iter().enumerate() {
        let atom2 = ring[(index + 1) % ring.len()];
        let Some(row) = topology
            .adjacency
            .neighbors_of(atom1)
            .iter()
            .find(|row| row.atom_index == atom2)
        else {
            continue;
        };
        let bond = &topology.bonds[row.bond.index()];
        if bond.order() != cosmolkit_model::BondOrder::Double
            || matches!(
                bond.stereo(),
                cosmolkit_model::BondStereo::None | cosmolkit_model::BondStereo::Any
            )
        {
            continue;
        }
        let Some([left_neighbor, right_neighbor]) = bond.stereo_atoms() else {
            continue;
        };
        let left_inside = ring.contains(&left_neighbor.index());
        let right_inside = ring.contains(&right_neighbor.index());
        let trans = match bond.stereo() {
            cosmolkit_model::BondStereo::Trans | cosmolkit_model::BondStereo::E => {
                left_inside == right_inside
            }
            _ => left_inside != right_inside,
        };
        if !trans {
            continue;
        }
        let left = coords[&ring[(index + ring.len() - 1) % ring.len()]];
        let reference = coords[&atom2];
        let interest = coords[&atom1];
        let dx = left[0] - reference[0];
        let dy = left[1] - reference[1];
        let denominator = dx * dx + dy * dy;
        let a = (dx * dx - dy * dy) / denominator;
        let b = 2.0 * dx * dy / denominator;
        coords.insert(
            atom1,
            [
                a * (interest[0] - reference[0]) + b * (interest[1] - reference[1]) + reference[0],
                b * (interest[0] - reference[0]) - a * (interest[1] - reference[1]) + reference[1],
            ],
        );
    }
}

#[derive(Debug, Clone, PartialEq)]
pub(crate) struct EmbeddedAtom {
    pub(crate) aid: usize,
    pub(crate) angle: f64,
    pub(crate) nbr1: Option<usize>,
    pub(crate) nbr2: Option<usize>,
    pub(crate) cis_trans_nbr: Option<usize>,
    pub(crate) ccw: bool,
    pub(crate) rot_dir: i32,
    pub(crate) loc: Point2,
    pub(crate) normal: Point2,
    pub(crate) neighs: Vec<usize>,
    pub(crate) density: f64,
    pub(crate) fixed: bool,
}

impl EmbeddedAtom {
    fn source_default() -> Self {
        // RDKit❗✔️: EmbeddedAtom() { neighs.clear(); }
        // RDKit❗✔️: unsigned int aid{0};
        // RDKit❗✔️: double angle{-1.0};
        // RDKit❗✔️: int nbr1{-1};
        // RDKit❗✔️: int nbr2{-1};
        // RDKit❗✔️: int CisTransNbr{-1};
        // RDKit❗✔️: bool ccw{true};
        // RDKit❗✔️: int rotDir{0};
        // RDKit❗✔️: double d_density{-1.0};
        // RDKit❗✔️: bool df_fixed{false};
        // RDKit❗✔️: double x{0.0};
        // RDKit❗✔️: double y{0.0};
        // Behavior: this is the value inserted by std::map::operator[] in
        // the pinned partial-fragment flip/permutation paths. The final two
        // lines are Point2D's in-class member initializers.
        // Complexity: one fixed-size value construction with one empty Vec,
        // matching the source default object plus empty neighbor vector.
        Self::at(0, [0.0, 0.0])
    }

    fn at(aid: usize, loc: Point2) -> Self {
        // RDKit❗✔️: EmbeddedAtom(unsigned int aid, const RDGeom::Point2D &pos)
        // RDKit❗✔️:     : aid(aid),
        // RDKit❗✔️:       angle(-1.0),
        // RDKit❗✔️:       nbr1(-1),
        // RDKit❗✔️:       nbr2(-1),
        // RDKit❗✔️:       CisTransNbr(-1),
        // RDKit❗✔️:       ccw(true),
        // RDKit❗✔️:       rotDir(0),
        // RDKit❗✔️:       d_density(-1.0),
        // RDKit❗✔️:       df_fixed(false) {
        // RDKit❗✔️:   loc = pos;
        // RDKit❗✔️: }
        Self {
            aid,
            angle: -1.0,
            nbr1: None,
            nbr2: None,
            cis_trans_nbr: None,
            ccw: true,
            rot_dir: 0,
            loc,
            normal: [0.0, 0.0],
            neighs: Vec::new(),
            density: -1.0,
            fixed: false,
        }
    }

    fn transform(&mut self, trans: Transform2D) {
        // RDKit❗✔️: void Transform(const RDGeom::Transform2D &trans) {
        // RDKit❗✔️:   RDGeom::Point2D temp = loc + normal;
        // RDKit❗✔️:   trans.TransformPoint(loc);
        // RDKit❗✔️:   trans.TransformPoint(temp);
        // RDKit❗✔️:   normal = temp - loc;
        // RDKit❗✔️: }
        let tip = [self.loc[0] + self.normal[0], self.loc[1] + self.normal[1]];
        self.loc = trans.transform_point(self.loc);
        let tip = trans.transform_point(tip);
        self.normal = [tip[0] - self.loc[0], tip[1] - self.loc[1]];
    }

    fn reflect(&mut self, a: Point2, b: Point2) {
        // RDKit❗✔️: void Reflect(const RDGeom::Point2D &loc1, const RDGeom::Point2D &loc2) {
        // RDKit❗✔️:   RDGeom::Point2D temp = loc + normal;
        // RDKit❗✔️:   loc = reflectPoint(loc, loc1, loc2);
        // RDKit❗✔️:   temp = reflectPoint(temp, loc1, loc2);
        // RDKit❗✔️:   normal = temp - loc;
        // RDKit❗✔️:   ccw = (!ccw);
        // RDKit❗✔️: }
        let tip = [self.loc[0] + self.normal[0], self.loc[1] + self.normal[1]];
        self.loc = reflect_point(self.loc, a, b);
        let tip = reflect_point(tip, a, b);
        self.normal = [tip[0] - self.loc[0], tip[1] - self.loc[1]];
        self.ccw = !self.ccw;
    }
}

#[derive(Debug)]
pub(crate) struct EmbeddedFrag<'a> {
    pub(crate) atoms: BTreeMap<usize, EmbeddedAtom>,
    pub(crate) attachment_points: Vec<usize>,
    pub(crate) done: bool,
    topology: &'a TopologyBlock,
    rings: &'a RingInfo,
}

pub(crate) fn orient_and_shift_fragments(
    fragments: &mut [EmbeddedFrag<'_>],
    canonicalize: bool,
    coordinate_map_size: Option<usize>,
) {
    // RDKit❗✔️:   if (!params.coordMap || !params.coordMap->size()) {
    // RDKit❗✔️:     if (params.canonOrient && efrags.size()) {
    // RDKit❗✔️:       // if we do not have any prespecified coordinates - canonicalize
    // RDKit❗✔️:       // the orientation of the fragment so that the longest axes fall
    // RDKit❗✔️:       // along the x-axis etc.
    // RDKit❗✔️:       for (auto &eri : efrags) {
    // RDKit❗✔️:         eri.canonicalizeOrientation();
    // RDKit❗✔️:       }
    // RDKit❗✔️:     }
    // RDKit❗✔️:   }
    // RDKit❗✔️:   DepictorLocal::_shiftCoords(efrags);
    // Behavior: an absent or present-but-empty coordinate map permits
    // canonical orientation; any nonempty map suppresses it. Fragment shifting
    // is unconditional after that decision.
    // Complexity: one linear canonicalization pass per fragment when enabled,
    // followed by one linear bounding-box/translation pass.
    if coordinate_map_size.is_none_or(|size| size == 0) && canonicalize && !fragments.is_empty() {
        for fragment in fragments.iter_mut() {
            fragment.canonicalize_orientation();
        }
    }
    shift_coordinates(fragments);
}

pub(crate) fn seed_coordinate_constraints<'a>(
    topology: &'a TopologyBlock,
    rings: &'a RingInfo,
    coordinate_map: Option<&PointMap>,
) -> Result<Option<EmbeddedFrag<'a>>, FragmentError> {
    // RDKit❗✔️:   // user-specified coordinates exist
    // RDKit❗✔️:   bool preSpec = false;
    // RDKit❗✔️:   // first embed any atoms for which the coordinates have been specified.
    // RDKit❗✔️:   if ((coordMap) && (coordMap->size() > 1)) {
    // RDKit❗✔️:     EmbeddedFrag efrag(&mol, *coordMap);
    // RDKit❗✔️:     // add this to the list of embedded fragments
    // RDKit❗✔️:     efrags.push_back(efrag);
    // RDKit❗✔️:     preSpec = true;
    // RDKit❗✔️:   }
    // Behavior: every supplied index is checked before any fragment is
    // published; only a map with more than one row becomes the prespecified
    // starting fragment. Absent, empty and singleton maps do not seed.
    // Complexity: one ordered validation traversal and, for a multi-row map,
    // the source-shaped ordered constructor traversal.
    let Some(coordinate_map) = coordinate_map else {
        return Ok(None);
    };
    validate_coordinate_constraints(topology, coordinate_map)?;
    (coordinate_map.len() > 1)
        .then(|| EmbeddedFrag::from_coord_map(topology, rings, coordinate_map))
        .transpose()
}

pub(crate) fn translate_single_coordinate_constraint(
    topology: &TopologyBlock,
    fragments: &mut [EmbeddedFrag<'_>],
    coordinate_map: Option<&PointMap>,
) -> Result<(), FragmentError> {
    // RDKit❗✔️:   // special case for a single-atom coordMap template
    // RDKit❗✔️:   if ((params.coordMap) && (params.coordMap->size() == 1)) {
    // RDKit❗✔️:     auto &conf = mol.getConformer(cid);
    // RDKit❗✔️:     auto cRef = params.coordMap->begin();
    // RDKit❗✔️:     const auto &confPos = conf.getAtomPos(cRef->first);
    // RDKit❗✔️:     auto refPos = cRef->second;
    // RDKit❗✔️:     refPos.x -= confPos.x;
    // RDKit❗✔️:     refPos.y -= confPos.y;
    // RDKit❗✔️:     for (auto i = 0u; i < conf.getNumAtoms(); ++i) {
    // RDKit❗✔️:       auto confPos = conf.getAtomPos(i);
    // RDKit❗✔️:       confPos.x += refPos.x;
    // RDKit❗✔️:       confPos.y += refPos.y;
    // RDKit❗✔️:       conf.setAtomPos(i, confPos);
    // RDKit❗✔️:     }
    // RDKit❗✔️:   }
    // Behavior: a singleton map translates every final fragment by one common
    // vector after layout and shifting. The last source-order write supplies
    // the mapped conformer position if a transient fragment overlap remains;
    // an unwritten conformer row has the source default origin.
    // Complexity: one fragment scan to reproduce copyCoordinate's last write,
    // then one allocation-free linear translation pass.
    let Some(coordinate_map) = coordinate_map else {
        return Ok(());
    };
    validate_coordinate_constraints(topology, coordinate_map)?;
    if coordinate_map.len() != 1 {
        return Ok(());
    }
    let (&anchor, &target) = coordinate_map.iter().next().expect("singleton map");
    let mut current = [0.0, 0.0];
    for fragment in fragments.iter() {
        if let Some(atom) = fragment.atoms.get(&anchor) {
            current = atom.loc;
        }
    }
    let shift = [target[0] - current[0], target[1] - current[1]];
    for fragment in fragments {
        fragment.translate(shift);
    }
    Ok(())
}

fn validate_coordinate_constraints(
    topology: &TopologyBlock,
    coordinate_map: &PointMap,
) -> Result<(), FragmentError> {
    // RDKit❗✔️:   unsigned int na = mol->getNumAtoms();
    // RDKit❗✔️:   for (const auto &cri : coordMap) {
    // RDKit❗✔️:     unsigned int aid = cri.first;
    // RDKit❗✔️:     CHECK_INVARIANT(aid < na, "");
    // RDKit❗✔️:     EmbeddedAtom eatom(aid, cri.second);
    // RDKit❗✔️:     eatom.neighs.clear();
    // RDKit❗✔️:     eatom.df_fixed = true;
    // RDKit❗✔️:     d_eatoms[aid] = eatom;
    // RDKit❗✔️:     d_done = false;
    // RDKit❗✔️:   }
    // Behavior: Rust reports the source invariant failure structurally and
    // performs no fragment mutation before the complete ordered map validates.
    // Complexity: one ordered map traversal with constant state.
    for &atom in coordinate_map.keys() {
        if atom >= topology.atoms.len() {
            return Err(FragmentError::AtomIndexOutOfRange {
                atom,
                atom_count: topology.atoms.len(),
            });
        }
    }
    Ok(())
}

fn shift_coordinates(fragments: &mut [EmbeddedFrag<'_>]) {
    // RDKit❗✔️: void _shiftCoords(std::list<EmbeddedFrag> &efrags) {
    // RDKit❗✔️:   // shift the coordinates if there are multiple fragments
    // RDKit❗✔️:   // so that the fragments do not overlap each other
    // RDKit❗✔️:   if (efrags.empty()) {
    // RDKit❗✔️:     return;
    // RDKit❗✔️:   }
    // RDKit❗✔️:   for (auto &efrag : efrags) {
    // RDKit❗✔️:     efrag.computeBox();
    // RDKit❗✔️:   }
    // RDKit❗✔️:   auto eri = efrags.begin();
    // RDKit❗✔️:   auto xmax = eri->getBoxPx();
    // RDKit❗✔️:   auto xmin = eri->getBoxNx();
    // RDKit❗✔️:   auto ymax = eri->getBoxPy();
    // RDKit❗✔️:   auto ymin = eri->getBoxNy();
    // RDKit❗✔️:
    // RDKit❗✔️:   ++eri;
    // RDKit❗✔️:   while (eri != efrags.end()) {
    // RDKit❗✔️:     bool xshift = true;
    // RDKit❗✔️:
    // RDKit❗✔️:     if (xmax + xmin > ymax + ymin) {
    // RDKit❗✔️:       xshift = false;
    // RDKit❗✔️:     }
    // RDKit❗✔️:     auto xn = eri->getBoxNx();
    // RDKit❗✔️:     auto xp = eri->getBoxPx();
    // RDKit❗✔️:     auto yn = eri->getBoxNy();
    // RDKit❗✔️:     auto yp = eri->getBoxPy();
    // RDKit❗✔️:     RDGeom::Point2D shift(0.0, 0.0);
    // RDKit❗✔️:     if (xshift) {
    // RDKit❗✔️:       shift.x = xmax + xn + 1.0;
    // RDKit❗✔️:       shift.y = 0.0;
    // RDKit❗✔️:       xmax += xp + xn + 1.0;
    // RDKit❗✔️:     } else {
    // RDKit❗✔️:       shift.x = 0.0;
    // RDKit❗✔️:       shift.y = ymax + yn + 1.0;
    // RDKit❗✔️:       ymax += yp + yn + 1.0;
    // RDKit❗✔️:     }
    // RDKit❗✔️:     eri->Translate(shift);
    // RDKit❗✔️:
    // RDKit❗✔️:     ++eri;
    // RDKit❗✔️:   }
    // RDKit❗✔️: }
    // Behavior: boxes and cumulative extents retain source list order and its
    // exact x-versus-y comparison; no packing heuristic is substituted.
    // Complexity: one box traversal per fragment and O(1) cross-fragment
    // state, matching the source asymptotic and avoiding an extra box vector.
    if fragments.is_empty() {
        return;
    }
    let [mut xmax, xmin, mut ymax, ymin] = fragments[0].compute_box();
    for fragment in &mut fragments[1..] {
        let [xp, xn, yp, yn] = fragment.compute_box();
        let xshift = xmax + xmin <= ymax + ymin;
        let shift = if xshift {
            let shift = [xmax + xn + 1.0, 0.0];
            xmax += xp + xn + 1.0;
            shift
        } else {
            let shift = [0.0, ymax + yn + 1.0];
            ymax += yp + yn + 1.0;
            shift
        };
        fragment.translate(shift);
    }
}

impl<'a> EmbeddedFrag<'a> {
    pub(crate) fn from_cis_trans_bond(
        bond_id: BondId,
        topology: &'a TopologyBlock,
        rings: &'a RingInfo,
    ) -> Result<Self, FragmentError> {
        // RDKit❗✔️: EmbeddedFrag::EmbeddedFrag(const RDKit::Bond *dblBond) {
        // RDKit❗✔️:   // Earlier embedding a cis/trans system meant to assign coordinates to the
        // RDKit❗✔️:   // atoms on the double bond as well as the neighboring atoms connected by the
        // RDKit❗✔️:   // single bond for which the cis/trans code has been specified. this causes
        // RDKit❗✔️:   // some ugliness in cases where these neighboring atoms are either part of a
        // RDKit❗✔️:   // different cis/trans system or a ring system. The function "merge" used to
        // RDKit❗✔️:   // deal with this ugliness. Now we will just embed the atoms on the double
        // RDKit❗✔️:   // bonds and mark at these atoms the direction in which the incoming single
        // RDKit❗✔️:   // bonds should go. Makes the merge function easier and address issue 171
        // RDKit❗✔️:   // simultaneously.
        // RDKit❗✔️:   PRECONDITION(dblBond, "");
        // RDKit❗✔️:   PRECONDITION(dblBond->getBondType() == RDKit::Bond::DOUBLE, "");
        // RDKit❗✔️:   auto stype = dblBond->getStereo();
        // RDKit❗✔️:   PRECONDITION(stype > RDKit::Bond::STEREOANY, "");
        // RDKit❗✔️:   const auto &nbrAtms = dblBond->getStereoAtoms();
        // RDKit❗✔️:   PRECONDITION(nbrAtms.size() == 2, "");
        // RDKit❗✔️:   dp_mol = &(dblBond->getOwningMol());
        // RDKit❗✔️:   auto begAtm = dblBond->getBeginAtomIdx();
        // RDKit❗✔️:   auto endAtm = dblBond->getEndAtomIdx();
        // RDKit❗✔️:   // the begin atom goes at the origin and the normal goes along -ve y-axis
        // RDKit❗✔️:   // to be rotate clock to add the cis/trans single bond
        // RDKit❗✔️:   EmbeddedAtom beatm;
        // RDKit❗✔️:   beatm.aid = begAtm;
        // RDKit❗✔️:   beatm.loc = RDGeom::Point2D(0.0, 0.0);
        // RDKit❗✔️:   beatm.nbr1 = endAtm;
        // RDKit❗✔️:   beatm.normal = RDGeom::Point2D(0.0, -1.0);
        // RDKit❗✔️:   beatm.ccw = false;
        // RDKit❗✔️:   beatm.CisTransNbr = nbrAtms[0];
        // RDKit❗✔️:   d_eatoms[begAtm] = beatm;
        // RDKit❗✔️:   // the end atom goes on the x-axis
        // RDKit❗✔️:   EmbeddedAtom eeatm;
        // RDKit❗✔️:   eeatm.aid = endAtm;
        // RDKit❗✔️:   eeatm.loc = RDGeom::Point2D(BOND_LEN, 0.0);
        // RDKit❗✔️:   eeatm.nbr1 = begAtm;
        // RDKit❗✔️:   eeatm.CisTransNbr = nbrAtms[1];
        // RDKit❗✔️:   if (stype == RDKit::Bond::STEREOZ || stype == RDKit::Bond::STEREOCIS) {
        // RDKit❗✔️:     eeatm.normal = RDGeom::Point2D(0.0, -1.0);
        // RDKit❗✔️:     eeatm.ccw = true;
        // RDKit❗✔️:   } else {
        // RDKit❗✔️:     eeatm.normal = RDGeom::Point2D(0.0, 1.0);
        // RDKit❗✔️:     eeatm.ccw = false;
        // RDKit❗✔️:   }
        // RDKit❗✔️:   d_eatoms[endAtm] = eeatm;
        // RDKit❗✔️:   d_done = false;
        // RDKit❗✔️: }
        // Behavior: only the two alkene endpoints are seeded. The designated
        // substituents remain pending neighbors, with source cis/trans normals.
        // Complexity: two embedded atom rows and no whole-topology clone.
        let bond =
            topology
                .bonds
                .get(bond_id.index())
                .ok_or(FragmentError::CisTransBondInvalid {
                    bond: bond_id.index(),
                })?;
        if bond.order() != BondOrder::Double
            || bond.stereo().rdkit_code() <= BondStereo::Any.rdkit_code()
        {
            return Err(FragmentError::CisTransBondInvalid {
                bond: bond_id.index(),
            });
        }
        let [left, right] = bond
            .stereo_atoms()
            .ok_or(FragmentError::CisTransBondInvalid {
                bond: bond_id.index(),
            })?;
        let begin = bond.begin().index();
        let end = bond.end().index();
        let mut first = EmbeddedAtom::at(begin, [0.0, 0.0]);
        first.nbr1 = Some(end);
        first.normal = [0.0, -1.0];
        first.ccw = false;
        first.cis_trans_nbr = Some(left.index());
        let mut second = EmbeddedAtom::at(end, [BOND_LEN, 0.0]);
        second.nbr1 = Some(begin);
        second.cis_trans_nbr = Some(right.index());
        if matches!(bond.stereo(), BondStereo::Z | BondStereo::Cis) {
            second.normal = [0.0, -1.0];
            second.ccw = true;
        } else {
            second.normal = [0.0, 1.0];
            second.ccw = false;
        }
        Ok(Self {
            atoms: BTreeMap::from([(begin, first), (end, second)]),
            attachment_points: Vec::new(),
            done: false,
            topology,
            rings,
        })
    }
    pub(crate) fn from_single(
        aid: usize,
        topology: &'a TopologyBlock,
        rings: &'a RingInfo,
    ) -> Result<Self, FragmentError> {
        // RDKit❗✔️: EmbeddedFrag::EmbeddedFrag(unsigned int aid, const RDKit::ROMol *mol) {
        // RDKit❗✔️:   PRECONDITION(mol, "");
        // RDKit❗✔️:   PRECONDITION(aid < mol->getNumAtoms(), "");
        // RDKit❗✔️:
        // RDKit❗✔️:   EmbeddedAtom eatm;
        // RDKit❗✔️:   eatm.aid = aid;
        // RDKit❗✔️:   RDGeom::Point2D org(0.0, 0.0);
        // RDKit❗✔️:   RDGeom::Point2D normal(1.0, 0.0);
        // RDKit❗✔️:   eatm.loc = org;
        // RDKit❗✔️:   eatm.normal = normal;
        // RDKit❗✔️:   eatm.angle = -1.0;
        // RDKit❗✔️:   eatm.ccw = true;
        // RDKit❗✔️:   eatm.neighs.clear();
        // RDKit❗✔️:   d_eatoms.clear();
        // RDKit❗✔️:   d_attachPts.clear();
        // RDKit❗✔️:   d_eatoms[aid] = eatm;
        // RDKit❗✔️:   d_done = false;
        // RDKit❗✔️:   dp_mol = mol;
        // RDKit❗✔️:   this->updateNewNeighs(aid);
        // RDKit❗✔️: }
        if aid >= topology.atoms.len() {
            return Err(FragmentError::AtomIndexOutOfRange {
                atom: aid,
                atom_count: topology.atoms.len(),
            });
        }
        let mut atom = EmbeddedAtom::at(aid, [0.0, 0.0]);
        atom.normal = [1.0, 0.0];
        let mut fragment = Self {
            atoms: BTreeMap::from([(aid, atom)]),
            attachment_points: Vec::new(),
            done: false,
            topology,
            rings,
        };
        fragment.update_new_neighbors(aid)?;
        Ok(fragment)
    }

    pub(crate) fn from_coord_map(
        topology: &'a TopologyBlock,
        rings: &'a RingInfo,
        coord_map: &PointMap,
    ) -> Result<Self, FragmentError> {
        // RDKit❗✔️: EmbeddedFrag::EmbeddedFrag(const RDKit::ROMol *mol,
        // RDKit❗✔️:                            const RDGeom::INT_POINT2D_MAP &coordMap) {
        // RDKit❗✔️:   PRECONDITION(mol, "");
        // RDKit❗✔️:   dp_mol = mol;
        // RDKit❗✔️:   d_eatoms.clear();
        // RDKit❗✔️:   d_attachPts.clear();
        // RDKit❗✔️:   unsigned int na = mol->getNumAtoms();
        // RDKit❗✔️:   for (const auto &cri : coordMap) {
        // RDKit❗✔️:     unsigned int aid = cri.first;
        // RDKit❗✔️:     CHECK_INVARIANT(aid < na, "");
        // RDKit❗✔️:     EmbeddedAtom eatom(aid, cri.second);
        // RDKit❗✔️:     eatom.neighs.clear();
        // RDKit❗✔️:     eatom.df_fixed = true;
        // RDKit❗✔️:     d_eatoms[aid] = eatom;
        // RDKit❗✔️:     d_done = false;
        // RDKit❗✔️:   }
        // RDKit❗✔️:   this->setupNewNeighs();
        // RDKit❗✔️:   this->setupAttachmentPoints();
        // RDKit❗✔️: }
        let mut atoms = BTreeMap::new();
        for (&aid, &loc) in coord_map {
            if aid >= topology.atoms.len() {
                return Err(FragmentError::AtomIndexOutOfRange {
                    atom: aid,
                    atom_count: topology.atoms.len(),
                });
            }
            let mut atom = EmbeddedAtom::at(aid, loc);
            atom.fixed = true;
            atoms.insert(aid, atom);
        }
        let mut fragment = Self {
            atoms,
            attachment_points: Vec::new(),
            done: false,
            topology,
            rings,
        };
        fragment.setup_new_neighbors()?;
        fragment.setup_attachment_points()?;
        Ok(fragment)
    }

    pub(crate) fn from_fused_rings(
        topology: &'a TopologyBlock,
        rings: &'a RingInfo,
        fused_rings: &[Vec<usize>],
        use_templates: bool,
        templates: &mut CoordinateTemplates,
    ) -> Result<Self, FragmentError> {
        // RDKit❗✔️: EmbeddedFrag::EmbeddedFrag(const RDKit::ROMol *mol,
        // RDKit❗✔️:                            const RDKit::VECT_INT_VECT &fusedRings,
        // RDKit❗✔️:                            bool useRingTemplates) {
        // RDKit❗✔️:   PRECONDITION(mol, "");
        // RDKit❗✔️:   dp_mol = mol;
        // RDKit❗✔️:   d_eatoms.clear();
        // RDKit❗✔️:   d_attachPts.clear();
        // RDKit❗✔️:   this->embedFusedRings(fusedRings, useRingTemplates);
        // RDKit❗✔️:   d_done = false;
        // RDKit❗✔️: }
        // Behavior: source constructor initializes a fresh fragment and then
        // embeds the complete ordered ring system before neighbor refresh.
        // Complexity: no topology copy; detached slices are borrowed.
        let mut fragment = Self {
            atoms: BTreeMap::new(),
            attachment_points: Vec::new(),
            done: false,
            topology,
            rings,
        };
        fragment.embed_fused_rings(fused_rings, use_templates, templates)?;
        Ok(fragment)
    }

    fn init_from_ring_coords(&mut self, ring: &[usize], coords: &PointMap) {
        // RDKit❗✔️: void EmbeddedFrag::initFromRingCoords(const RDKit::INT_VECT &ring,
        // RDKit❗✔️:                                       const RDGeom::INT_POINT2D_MAP &nringMap) {
        // RDKit❗✔️:   double largestAngle = M_PI * (1 - (2.0 / ring.size()));
        // RDKit❗✔️:   auto prev = ring.back();
        // RDKit❗✔️:   unsigned int cnt = 0;
        // RDKit❗✔️:   for (auto ai : ring) {
        // RDKit❗✔️:     EmbeddedAtom eatm;
        // RDKit❗✔️:     eatm.loc = nringMap.at(ai);
        // RDKit❗✔️:     eatm.aid = ai;
        // RDKit❗✔️:     eatm.angle = largestAngle;
        // RDKit❗✔️:     eatm.nbr1 = prev;
        // RDKit❗✔️:     if (cnt) {
        // RDKit❗✔️:       d_eatoms[prev].nbr2 = ai;
        // RDKit❗✔️:     }
        // RDKit❗✔️:     d_eatoms[ai] = eatm;
        // RDKit❗✔️:     prev = ai;
        // RDKit❗✔️:     cnt++;
        // RDKit❗✔️:   }
        // RDKit❗✔️:   d_eatoms[prev].nbr2 = ring.front();
        // RDKit❗✔️: }
        // Behavior: ordered cyclic predecessor/successor and source angle.
        // Complexity: one map insertion per ring atom.
        let angle = PI * (1.0 - 2.0 / ring.len() as f64);
        let mut previous = *ring.last().expect("ring decomposition has nonempty rings");
        for (index, &atom_id) in ring.iter().enumerate() {
            let mut atom = EmbeddedAtom::at(atom_id, coords[&atom_id]);
            atom.angle = angle;
            atom.nbr1 = Some(previous);
            if index != 0 {
                self.atoms
                    .get_mut(&previous)
                    .expect("previous ring atom")
                    .nbr2 = Some(atom_id);
            }
            self.atoms.insert(atom_id, atom);
            previous = atom_id;
        }
        self.atoms.get_mut(&previous).expect("last ring atom").nbr2 = Some(ring[0]);
    }

    fn merge_ring(&mut self, other: &Self, common_count: usize, pin_atoms: &[usize]) {
        // RDKit❗✔️: void EmbeddedFrag::mergeRing(const EmbeddedFrag &embRing, unsigned int nCommon,
        // RDKit❗✔️:                              const RDKit::INT_VECT &pinAtoms) {
        // RDKit❗✔️:   const auto &oatoms = embRing.GetEmbeddedAtoms();
        // RDKit❗✔️:   for (const auto &ori : oatoms) {
        // RDKit❗✔️:     auto aid = ori.first;
        // RDKit❗✔️:     if (d_eatoms.find(aid) == d_eatoms.end()) {
        // RDKit❗✔️:       d_eatoms[aid] = ori.second;
        // RDKit❗✔️:     } else {
        // RDKit❗✔️:       if (nCommon <= 2) {
        // RDKit❗✔️:         if (std::find(pinAtoms.begin(), pinAtoms.end(), aid) !=
        // RDKit❗✔️:             pinAtoms.end()) {
        // RDKit❗✔️:           d_eatoms[aid].angle += ori.second.angle;
        // RDKit❗✔️:           if (d_eatoms[aid].nbr1 == ori.second.nbr1) {
        // RDKit❗✔️:             d_eatoms[aid].nbr1 = ori.second.nbr2;
        // RDKit❗✔️:           } else if (d_eatoms[aid].nbr1 == ori.second.nbr2) {
        // RDKit❗✔️:             d_eatoms[aid].nbr1 = ori.second.nbr1;
        // RDKit❗✔️:           } else if (d_eatoms[aid].nbr2 == ori.second.nbr1) {
        // RDKit❗✔️:             d_eatoms[aid].nbr2 = ori.second.nbr2;
        // RDKit❗✔️:           } else if (d_eatoms[aid].nbr2 == ori.second.nbr2) {
        // RDKit❗✔️:             d_eatoms[aid].nbr2 = ori.second.nbr1;
        // RDKit❗✔️:           }
        // RDKit❗✔️:         }
        // RDKit❗✔️:       }
        // RDKit❗✔️:     }
        // RDKit❗✔️:   }
        // RDKit❗✔️: }
        // Behavior: bridged overlaps (>2) retain existing neighbor state;
        // only pinned one-/two-atom overlaps accumulate angles.
        // Complexity: BTreeMap membership is logarithmic rather than source
        // ordered-map logarithmic, and each atom is visited once.
        for (&id, incoming) in &other.atoms {
            if let Some(current) = self.atoms.get_mut(&id) {
                if common_count <= 2 && pin_atoms.contains(&id) {
                    current.angle += incoming.angle;
                    if current.nbr1 == incoming.nbr1 {
                        current.nbr1 = incoming.nbr2;
                    } else if current.nbr1 == incoming.nbr2 {
                        current.nbr1 = incoming.nbr1;
                    } else if current.nbr2 == incoming.nbr1 {
                        current.nbr2 = incoming.nbr2;
                    } else if current.nbr2 == incoming.nbr2 {
                        current.nbr2 = incoming.nbr1;
                    }
                }
            } else {
                self.atoms.insert(id, incoming.clone());
            }
        }
    }

    fn embed_fused_rings(
        &mut self,
        fused_rings: &[Vec<usize>],
        use_templates: bool,
        templates: &mut CoordinateTemplates,
    ) -> Result<(), FragmentError> {
        // RDKit❗❗: void EmbeddedFrag::embedFusedRings(const RDKit::VECT_INT_VECT &fusedRings,
        // RDKit❗❗:                                    bool useRingTemplates) {
        // RDKit❗❗:   PRECONDITION(dp_mol, "");
        // RDKit❗❗:   RDKit::INT_VECT funion;
        // RDKit❗❗:   if (useRingTemplates &&
        // RDKit❗❗:       (fusedRings.size() > 1 ||
        // RDKit❗❗:        (fusedRings.size() == 1 && fusedRings[0].size() > 8))) {
        // RDKit❗❗:     RDKit::Union(fusedRings, funion);
        // RDKit❗❗:     bool found_template = matchToTemplate(funion, fusedRings.size());
        // RDKit❗❗:     if (found_template) {
        // RDKit❗❗:       return;
        // RDKit❗❗:     }
        // RDKit❗❗:   }
        // RDKit❗❗:   std::vector<RDGeom::INT_POINT2D_MAP> coords;
        // RDKit❗❗:   coords.reserve(fusedRings.size());
        // RDKit❗❗:   for (const auto &ring : fusedRings) {
        // RDKit❗❗:     auto ring_coords = embedRing(ring);
        // RDKit❗❗:     mirrorTransRingAtoms(*dp_mol, ring, ring_coords);
        // RDKit❗❗:     coords.push_back(ring_coords);
        // RDKit❗❗:   }
        // RDKit❗❗:   RDKit::INT_VECT doneRings;
        // RDKit❗❗:   if (useRingTemplates) {
        // RDKit❗❗:     RDKit::INT_VECT coreRingsIds;
        // RDKit❗❗:     auto coreRings = findCoreRings(fusedRings, coreRingsIds, *dp_mol);
        // RDKit❗❗:     if (coreRings.size() > 1 && coreRings.size() < fusedRings.size()) {
        // RDKit❗❗:       RDKit::Union(coreRings, funion);
        // RDKit❗❗:       bool found_template = matchToTemplate(funion, coreRings.size());
        // RDKit❗❗:       if (found_template) {
        // RDKit❗❗:         doneRings = coreRingsIds;
        // RDKit❗❗:       }
        // RDKit❗❗:     }
        // RDKit❗❗:   }
        // RDKit❗❗:   if (doneRings.empty()) {
        // RDKit❗❗:     auto firstRingId = pickFirstRingToEmbed(*dp_mol, fusedRings);
        // RDKit❗❗:     this->initFromRingCoords(fusedRings[firstRingId], coords[firstRingId]);
        // RDKit❗❗:     doneRings.push_back(firstRingId);
        // RDKit❗❗:   }
        // RDKit❗❗:   RDKit::Union(fusedRings, funion);
        // RDKit❗❗:   while (d_eatoms.size() < funion.size()) {
        // RDKit❗❗:     int nextId;
        // RDKit❗❗:     auto commonAtomIds = findNextRingToEmbed(doneRings, fusedRings, nextId);
        // RDKit❗❗:     RDGeom::Transform2D trans;
        // RDKit❗❗:     EmbeddedFrag embRing;
        // RDKit❗❗:     embRing.initFromRingCoords(fusedRings[nextId], coords[nextId]);
        // RDKit❗❗:     RDKit::INT_VECT pinAtoms;
        // RDKit❗❗:     if (commonAtomIds.size() == 1) {
        // RDKit❗❗:       trans.assign(this->computeOneAtomTrans(commonAtomIds[0], embRing));
        // RDKit❗❗:       embRing.Transform(trans);
        // RDKit❗❗:       pinAtoms.push_back(commonAtomIds.front());
        // RDKit❗❗:     } else {
        // RDKit❗❗:       auto aid1 = commonAtomIds.front();
        // RDKit❗❗:       auto aid2 = commonAtomIds.back();
        // RDKit❗❗:       pinAtoms.push_back(aid1);
        // RDKit❗❗:       pinAtoms.push_back(aid2);
        // RDKit❗❗:       trans.assign(this->computeTwoAtomTrans(aid1, aid2, coords[nextId]));
        // RDKit❗❗:       embRing.Transform(trans);
        // RDKit❗❗:       reflectIfNecessaryDensity(embRing, aid1, aid2);
        // RDKit❗❗:     }
        // RDKit❗❗:     this->mergeRing(embRing, commonAtomIds.size(), pinAtoms);
        // RDKit❗❗:     doneRings.push_back(nextId);
        // RDKit❗❗:   }
        // RDKit❗❗: }
        // Behavior: source-order full/core template attempts, seed choice,
        // cyclic mirror, overlap preference, transform/reflection and merge.
        // Complexity: one geometry map per ring and repeated overlap scans as
        // in source; BTreeMap versus std::map has equivalent logarithmic cost.
        let union = ring_union(fused_rings);
        if use_templates
            && (fused_rings.len() > 1 || fused_rings.first().is_some_and(|ring| ring.len() > 8))
            && self.match_to_template(&union, fused_rings.len(), templates)?
        {
            return Ok(());
        }
        let coordinates: Vec<_> = fused_rings
            .iter()
            .map(|ring| {
                let mut points = embed_ring(ring);
                mirror_trans_ring_atoms(self.topology, ring, &mut points);
                points
            })
            .collect();
        let mut done = Vec::new();
        if use_templates {
            let core_ids = find_core_rings(self.topology, fused_rings);
            if core_ids.len() > 1 && core_ids.len() < fused_rings.len() {
                let core_rings: Vec<_> = core_ids
                    .iter()
                    .map(|&index| fused_rings[index].clone())
                    .collect();
                if self.match_to_template(&ring_union(&core_rings), core_ids.len(), templates)? {
                    done = core_ids;
                }
            }
        }
        if done.is_empty() {
            let first = pick_first_ring_to_embed(self.topology, fused_rings);
            self.init_from_ring_coords(&fused_rings[first], &coordinates[first]);
            done.push(first);
        }
        while self.atoms.len() < union.len() {
            let (next, common) = find_next_ring_to_embed(&done, fused_rings);
            let mut incoming = Self {
                atoms: BTreeMap::new(),
                attachment_points: Vec::new(),
                done: false,
                topology: self.topology,
                rings: self.rings,
            };
            incoming.init_from_ring_coords(&fused_rings[next], &coordinates[next]);
            let pins = if common.len() == 1 {
                let transform = self.compute_one_atom_trans(common[0], &incoming)?;
                incoming.transform(transform);
                vec![common[0]]
            } else {
                let first = common[0];
                let last = *common.last().expect("connected ring overlap");
                let transform = self.compute_two_atom_trans(first, last, &incoming)?;
                incoming.transform(transform);
                self.reflect_if_necessary_density(&mut incoming, first, last)?;
                vec![first, last]
            };
            self.merge_ring(&incoming, common.len(), &pins);
            done.push(next);
        }
        Ok(())
    }

    pub(crate) fn match_to_template(
        &mut self,
        ring_system_atoms: &[usize],
        ring_count: usize,
        templates: &mut CoordinateTemplates,
    ) -> Result<bool, FragmentError> {
        // BEGIN RDKIT CPP FUNCTION EmbeddedFrag::matchToTemplate
        // RDKit❗❗: bool EmbeddedFrag::matchToTemplate(const RDKit::INT_VECT &ringSystemAtoms,
        // RDKit❗❗:                                    unsigned int ring_count) {
        // RDKit❗❗:   CoordinateTemplates &coordinate_templates =
        // RDKit❗❗:       CoordinateTemplates::getRingSystemTemplates();
        // RDKit❗❗:
        // RDKit❗❗:   // only look for an exact match to the ring system because our method of
        // RDKit❗❗:   // completing rings from a template isn't reliably better than not using
        // RDKit❗❗:   // a template at all
        // RDKit❗❗:   if (!coordinate_templates.hasTemplateOfSize(ringSystemAtoms.size())) {
        // RDKit❗❗:     return false;
        // RDKit❗❗:   }
        // RDKit❗❗:
        // RDKit❗❗:   // make a mol out of the induced subgraph using the ring system atoms
        // RDKit❗❗:   RDKit::RWMol rs_mol(*dp_mol, true);
        // RDKit❗❗:
        // RDKit❗❗:   boost::dynamic_bitset<> rs_atoms(dp_mol->getNumAtoms());
        // RDKit❗❗:   for (auto aidx : ringSystemAtoms) {
        // RDKit❗❗:     rs_atoms.set(aidx);
        // RDKit❗❗:   }
        // RDKit❗❗:
        // RDKit❗❗:   constexpr int DUMMY_ATOMIC_NUM = 200;
        // RDKit❗❗:   for (auto &at : rs_mol.atoms()) {
        // RDKit❗❗:     if (!rs_atoms.test(at->getIdx())) {
        // RDKit❗❗:       at->setAtomicNum(DUMMY_ATOMIC_NUM);
        // RDKit❗❗:     }
        // RDKit❗❗:   }
        // RDKit❗❗:   auto numBonds = rs_mol.getNumBonds();
        // RDKit❗❗:   for (auto bnd : rs_mol.bonds()) {
        // RDKit❗❗:     if (!rs_atoms.test(bnd->getBeginAtomIdx()) ||
        // RDKit❗❗:         !rs_atoms.test(bnd->getEndAtomIdx())) {
        // RDKit❗❗:       --numBonds;
        // RDKit❗❗:     }
        // RDKit❗❗:   }
        // RDKit❗❗:
        // RDKit❗❗:   // find template that this mol matches to, if any
        // RDKit❗❗:   RDKit::MatchVectType match;
        // RDKit❗❗:   std::shared_ptr<RDKit::ROMol> template_mol(nullptr);
        // RDKit❗❗:   for (const auto &mol :
        // RDKit❗❗:        coordinate_templates.getMatchingTemplates(ringSystemAtoms.size())) {
        // RDKit❗❗:     // To reduce how often we have to do substructure matches, check ring info
        // RDKit❗❗:     // and bond count first
        // RDKit❗❗:     if (mol->getNumBonds() != numBonds) {
        // RDKit❗❗:       continue;
        // RDKit❗❗:     } else if (mol->getRingInfo()->numRings() != ring_count) {
        // RDKit❗❗:       continue;
        // RDKit❗❗:     }
        // RDKit❗❗:     // also check if the mol atoms have the same connectivity as the template
        // RDKit❗❗: #ifdef _MSC_VER
        // RDKit❗❗:     // MSVC++ doesn't like implicitly capturing constexpr variables, this is a
        // RDKit❗❗:     // bug
        // RDKit❗❗:     auto degreeCounts = [DUMMY_ATOMIC_NUM](const RDKit::ROMol &mol) {
        // RDKit❗❗: #else
        // RDKit❗❗:     // clang generates warnings if you explicitly capture a constexpr variable
        // RDKit❗❗:     auto degreeCounts = [](const RDKit::ROMol &mol) {
        // RDKit❗❗: #endif
        // RDKit❗❗:       std::array<int, 5> degrees_count({0, 0, 0, 0, 0});
        // RDKit❗❗:       for (auto atom : mol.atoms()) {
        // RDKit❗❗:         if (atom->getAtomicNum() == DUMMY_ATOMIC_NUM) {
        // RDKit❗❗:           continue;
        // RDKit❗❗:         }
        // RDKit❗❗:         auto degree = 0u;
        // RDKit❗❗:         for (auto nbr : mol.atomNeighbors(atom)) {
        // RDKit❗❗:           if (nbr->getAtomicNum() != DUMMY_ATOMIC_NUM) {
        // RDKit❗❗:             ++degree;
        // RDKit❗❗:             if (degree == 4) {
        // RDKit❗❗:               break;
        // RDKit❗❗:             }
        // RDKit❗❗:           }
        // RDKit❗❗:         }
        // RDKit❗❗:         degrees_count[degree]++;
        // RDKit❗❗:       }
        // RDKit❗❗:       return degrees_count;
        // RDKit❗❗:     };
        // RDKit❗❗:     if (degreeCounts(rs_mol) != degreeCounts(*mol)) {
        // RDKit❗❗:       continue;
        // RDKit❗❗:     }
        // RDKit❗❗:     RDKit::SubstructMatchParameters params;
        // RDKit❗❗:     params.maxMatches = 1;
        // RDKit❗❗:     auto matches = RDKit::SubstructMatch(rs_mol, *mol, params);
        // RDKit❗❗:     if (!matches.empty()) {
        // RDKit❗❗:       if (checkStereoChemistry(rs_mol, *mol, matches[0])) {
        // RDKit❗❗:         match = matches[0];
        // RDKit❗❗:         template_mol = mol;
        // RDKit❗❗:         break;
        // RDKit❗❗:       }
        // RDKit❗❗:     }
        // RDKit❗❗:   }
        // RDKit❗❗:   if (!template_mol) {
        // RDKit❗❗:     return false;
        // RDKit❗❗:   }
        // RDKit❗❗:
        // RDKit❗❗:   // copy over new coordinates
        // RDKit❗❗:   const auto &conf = template_mol->getConformer();
        // RDKit❗❗:   for (auto &[template_aidx, rs_aidx] : match) {
        // RDKit❗❗:     EmbeddedAtom new_at(rs_aidx, conf.getAtomPos(template_aidx));
        // RDKit❗❗:     new_at.df_fixed = true;
        // RDKit❗❗:     d_eatoms.emplace(rs_aidx, new_at);
        // RDKit❗❗:   }
        // RDKit❗❗:   this->setupNewNeighs();
        // RDKit❗❗:   this->setupAttachmentPoints();
        // RDKit❗❗:   return true;
        // RDKit❗❗: }
        // END RDKIT CPP FUNCTION EmbeddedFrag::matchToTemplate
        // Behavior: source-order first acceptable candidate; exact parity is
        // still pending the owning template-match regression gate.
        // Complexity: the target keeps one O(V) sentinel vector instead of
        // cloning a mol, but matcher context and candidate costs require the
        // full-gate review before a performance equivalence claim.
        if !templates.has_template_of_size(ring_system_atoms.len()) {
            return Ok(false);
        }

        let mut in_ring_system = vec![false; self.topology.atoms.len()];
        for &index in ring_system_atoms {
            if index >= in_ring_system.len() {
                return Err(FragmentError::AtomIndexOutOfRange {
                    atom: index,
                    atom_count: in_ring_system.len(),
                });
            }
            in_ring_system[index] = true;
        }
        // RDKit sets non-ring atoms to atomic number 200 on a temporary mol.
        // The search-only projection keeps canonical Element values valid.
        let atomic_numbers: Vec<_> = in_ring_system
            .iter()
            .map(|&included| (!included).then_some(200))
            .collect();
        let num_bonds = self
            .topology
            .bonds
            .iter()
            .filter(|bond| {
                in_ring_system[bond.begin().index()] && in_ring_system[bond.end().index()]
            })
            .count();
        let target_degrees = ring_system_degree_counts(self.topology, &in_ring_system);
        let coordinates = cosmolkit_model::CoordinateBlock::default();
        let target = SearchTarget::new(
            self.topology,
            &coordinates,
            &self.topology.stereo_groups,
            None,
            None,
        )
        .with_atomic_number_overrides(&atomic_numbers);

        for template in templates.matching_templates(ring_system_atoms.len()) {
            if template.query.num_bonds() != num_bonds
                || template.rings.num_rings() != ring_count
                || template_degree_counts(template) != target_degrees
            {
                continue;
            }
            let matches = try_get_substruct_matches_with_params(
                &target,
                &template.query,
                &SubstructMatchParams {
                    max_matches: 1,
                    ..Default::default()
                },
            )?;
            let Some(first) = matches.first() else {
                continue;
            };
            if !check_template_stereo(self.topology, template, &first.atom_mapping) {
                continue;
            }
            for (template_index, &target_index) in first.atom_mapping.iter().enumerate() {
                let mut atom =
                    EmbeddedAtom::at(target_index, template_point(template, template_index));
                atom.fixed = true;
                self.atoms.entry(target_index).or_insert(atom);
            }
            self.setup_new_neighbors()?;
            self.setup_attachment_points()?;
            return Ok(true);
        }
        Ok(false)
    }

    pub(crate) fn update_new_neighbors(&mut self, aid: usize) -> Result<(), FragmentError> {
        // RDKit❗✔️: void EmbeddedFrag::updateNewNeighs(
        // RDKit❗✔️:     unsigned int aid) {  //, const RDKit::ROMol *mol) {
        // RDKit❗✔️:   PRECONDITION(dp_mol, "");
        // RDKit❗✔️:
        // RDKit❗✔️:   d_eatoms[aid].neighs.clear();
        // RDKit❗✔️:   RDKit::INT_VECT hIndices;
        // RDKit❗✔️:   for (const auto nbr : dp_mol->atomNeighbors(dp_mol->getAtomWithIdx(aid))) {
        // RDKit❗✔️:     if (d_eatoms.find(nbr->getIdx()) == d_eatoms.end()) {
        // RDKit❗✔️:       if (dp_mol->getAtomWithIdx(nbr->getIdx())->getAtomicNum() != 1) {
        // RDKit❗✔️:         d_eatoms[aid].neighs.push_back(nbr->getIdx());
        // RDKit❗✔️:       } else {
        // RDKit❗✔️:         hIndices.push_back(nbr->getIdx());
        // RDKit❗✔️:       }
        // RDKit❗✔️:     }
        // RDKit❗✔️:   }
        // RDKit❗✔️:   d_eatoms[aid].neighs.insert(d_eatoms[aid].neighs.end(), hIndices.begin(),
        // RDKit❗✔️:                               hIndices.end());
        // RDKit❗✔️:
        // RDKit❗✔️:   auto deg = getDepictDegree(dp_mol->getAtomWithIdx(aid));
        // RDKit❗✔️:   if ((d_eatoms[aid].neighs.size() > 0) &&
        // RDKit❗✔️:       ((deg < 4) || (d_eatoms[aid].neighs.size() < 3))) {
        // RDKit❗✔️:     d_eatoms[aid].neighs = rankAtomsByRank(*dp_mol, d_eatoms[aid].neighs);
        // RDKit❗✔️:   } else if ((deg >= 4) && (d_eatoms[aid].neighs.size() >= 3)) {
        // RDKit❗✔️:     d_eatoms[aid].neighs = setNbrOrder(aid, d_eatoms[aid].neighs, *dp_mol);
        // RDKit❗✔️:   }
        // RDKit❗✔️:
        // RDKit❗✔️:   if (d_eatoms[aid].neighs.size() > 0) {
        // RDKit❗✔️:     if (std::find(d_attachPts.begin(), d_attachPts.end(),
        // RDKit❗✔️:                   static_cast<int>(aid)) == d_attachPts.end()) {
        // RDKit❗✔️:       d_attachPts.push_back(aid);
        // RDKit❗✔️:     }
        // RDKit❗✔️:   }
        // RDKit❗✔️: }
        if aid >= self.topology.atoms.len() {
            return Err(FragmentError::AtomIndexOutOfRange {
                atom: aid,
                atom_count: self.topology.atoms.len(),
            });
        }
        if !self.atoms.contains_key(&aid) {
            return Err(FragmentError::AtomNotEmbedded { atom: aid });
        }
        let mut neighbors = Vec::new();
        let mut hydrogens = Vec::new();
        for neighbor in self.topology.adjacency.neighbors_of(aid) {
            let index = neighbor.atom_index;
            if !self.atoms.contains_key(&index) {
                if self.topology.atoms[index].atomic_number() == 1 {
                    hydrogens.push(index);
                } else {
                    neighbors.push(index);
                }
            }
        }
        neighbors.extend(hydrogens);
        let degree = self.topology.adjacency.neighbors_of(aid).len();
        if !neighbors.is_empty() && (degree < 4 || neighbors.len() < 3) {
            neighbors = rank_atoms_by_rank(self.topology, &neighbors, true)?;
        } else if degree >= 4 && neighbors.len() >= 3 {
            neighbors = set_neighbor_order(self.topology, aid, &neighbors)?;
        }
        self.atoms.get_mut(&aid).expect("checked above").neighs = neighbors;
        if !self.atoms[&aid].neighs.is_empty() && !self.attachment_points.contains(&aid) {
            self.attachment_points.push(aid);
        }
        Ok(())
    }

    pub(crate) fn setup_new_neighbors(&mut self) -> Result<(), FragmentError> {
        // RDKit❗✔️: void EmbeddedFrag::setupNewNeighs() {  // const RDKit::ROMol *mol) {
        // RDKit❗✔️:   PRECONDITION(dp_mol, "");
        // RDKit❗✔️:
        // RDKit❗✔️:   d_attachPts.clear();
        // RDKit❗✔️:   for (const auto &eci : d_eatoms) {
        // RDKit❗✔️:     this->updateNewNeighs(eci.first);
        // RDKit❗✔️:   }
        // RDKit❗✔️:   d_attachPts = rankAtomsByRank(*dp_mol, d_attachPts);
        // RDKit❗✔️: }
        self.attachment_points.clear();
        let keys: Vec<_> = self.atoms.keys().copied().collect();
        for aid in keys {
            self.update_new_neighbors(aid)?;
        }
        self.attachment_points = rank_atoms_by_rank(self.topology, &self.attachment_points, true)?;
        Ok(())
    }

    pub(crate) fn find_neighbor(&self, aid: usize) -> Result<Option<usize>, FragmentError> {
        // RDKit❗✔️: int EmbeddedFrag::findNeighbor(
        // RDKit❗✔️:     unsigned int aid) {  //, const RDKit::ROMol *mol) {
        // RDKit❗✔️:   PRECONDITION(dp_mol, "");
        // RDKit❗✔️:
        // RDKit❗✔️:   for (const auto nbr : dp_mol->atomNeighbors(dp_mol->getAtomWithIdx(aid))) {
        // RDKit❗✔️:     if (d_eatoms.find(nbr->getIdx()) != d_eatoms.end()) {
        // RDKit❗✔️:       return nbr->getIdx();
        // RDKit❗✔️:     }
        // RDKit❗✔️:   }
        // RDKit❗✔️:   return -1;
        // RDKit❗✔️: }
        if aid >= self.topology.atoms.len() {
            return Err(FragmentError::AtomIndexOutOfRange {
                atom: aid,
                atom_count: self.topology.atoms.len(),
            });
        }
        Ok(self
            .topology
            .adjacency
            .neighbors_of(aid)
            .iter()
            .find_map(|neighbor| {
                self.atoms
                    .contains_key(&neighbor.atom_index)
                    .then_some(neighbor.atom_index)
            }))
    }

    pub(crate) fn setup_attachment_points(&mut self) -> Result<(), FragmentError> {
        // RDKit❗✔️: void EmbeddedFrag::setupAttachmentPoints() {
        // RDKit❗✔️:   for (auto dai : d_attachPts) {
        // RDKit❗✔️:     RDKit::INT_VECT doneNbrs;
        // RDKit❗✔️:     const auto &enbrs = d_eatoms[dai].neighs;
        // RDKit❗✔️:     for (const auto nbrAtom :
        // RDKit❗✔️:          dp_mol->atomNeighbors(dp_mol->getAtomWithIdx(dai))) {
        // RDKit❗✔️:       if (std::find(enbrs.begin(), enbrs.end(),
        // RDKit❗✔️:                     static_cast<int>(nbrAtom->getIdx())) == enbrs.end()) {
        // RDKit❗✔️:         doneNbrs.push_back(nbrAtom->getIdx());
        // RDKit❗✔️:       }
        // RDKit❗✔️:     }
        // RDKit❗✔️:     if (doneNbrs.empty()) {
        // RDKit❗✔️:       d_eatoms[dai].normal = RDGeom::Point2D(1., 0.);
        // RDKit❗✔️:       d_eatoms[dai].angle = -1.;
        // RDKit❗✔️:     } else if (doneNbrs.size() == 1) {
        // RDKit❗✔️:       auto nbid = doneNbrs.front();
        // RDKit❗✔️:       d_eatoms[dai].nbr1 = nbid;
        // RDKit❗✔️:       d_eatoms[dai].normal =
        // RDKit❗✔️:           computeNormal(d_eatoms[dai].loc, d_eatoms[nbid].loc);
        // RDKit❗✔️:     } else if (doneNbrs.size() == 2) {
        // RDKit❗✔️:       auto nb1 = doneNbrs[0];
        // RDKit❗✔️:       auto nb2 = doneNbrs[1];
        // RDKit❗✔️:       d_eatoms[dai].nbr1 = nb1;
        // RDKit❗✔️:       d_eatoms[dai].nbr2 = nb2;
        // RDKit❗✔️:       d_eatoms[dai].angle =
        // RDKit❗✔️:           computeAngle(d_eatoms[dai].loc, d_eatoms[nb1].loc, d_eatoms[nb2].loc);
        // RDKit❗✔️:     } else if (doneNbrs.size() >= 3) {
        // RDKit❗✔️:       this->computeNbrsAndAng(dai, doneNbrs);
        // RDKit❗✔️:     }
        // RDKit❗✔️:   }
        // RDKit❗✔️: }
        let attachment_points = self.attachment_points.clone();
        for aid in attachment_points {
            let atom = self
                .atoms
                .get(&aid)
                .ok_or(FragmentError::AtomNotEmbedded { atom: aid })?;
            let done_neighbors: Vec<_> = self
                .topology
                .adjacency
                .neighbors_of(aid)
                .iter()
                .map(|neighbor| neighbor.atom_index)
                .filter(|neighbor| !atom.neighs.contains(neighbor))
                .collect();
            match done_neighbors.as_slice() {
                [] => {
                    let atom = self.atoms.get_mut(&aid).expect("checked above");
                    atom.normal = [1.0, 0.0];
                    atom.angle = -1.0;
                }
                [first] => {
                    let center = self.atoms[&aid].loc;
                    let other = self
                        .atoms
                        .get(first)
                        .ok_or(FragmentError::AtomNotEmbedded { atom: *first })?
                        .loc;
                    let atom = self.atoms.get_mut(&aid).expect("checked above");
                    atom.nbr1 = Some(*first);
                    atom.normal = compute_normal(center, other)?;
                }
                [first, second] => {
                    let center = self.atoms[&aid].loc;
                    let a = self
                        .atoms
                        .get(first)
                        .ok_or(FragmentError::AtomNotEmbedded { atom: *first })?
                        .loc;
                    let b = self
                        .atoms
                        .get(second)
                        .ok_or(FragmentError::AtomNotEmbedded { atom: *second })?
                        .loc;
                    let angle = compute_angle(center, a, b)?;
                    let atom = self.atoms.get_mut(&aid).expect("checked above");
                    atom.nbr1 = Some(*first);
                    atom.nbr2 = Some(*second);
                    atom.angle = angle;
                }
                _ => self.compute_nbrs_and_ang(aid, &done_neighbors)?,
            }
        }
        Ok(())
    }

    fn compute_nbrs_and_ang(
        &mut self,
        aid: usize,
        done_neighbors: &[usize],
    ) -> Result<(), FragmentError> {
        // RDKit❗✔️: void EmbeddedFrag::computeNbrsAndAng(unsigned int aid,
        // RDKit❗✔️:                                      const RDKit::INT_VECT &doneNbrs) {
        // RDKit❗✔️:   PRECONDITION(dp_mol, "");
        // RDKit❗✔️:   PRECONDITION(aid < dp_mol->getNumAtoms(), "");
        // RDKit❗✔️:   PRECONDITION(doneNbrs.size() >= 3, "");
        // RDKit❗✔️:   std::list<DOUBLE_INT_PAIR> anglePairs;
        // RDKit❗✔️:   double ang;
        // RDKit❗✔️:   for (auto nbi1 = doneNbrs.begin(); nbi1 != doneNbrs.end(); ++nbi1) {
        // RDKit❗✔️:     auto nbi3 = nbi1;
        // RDKit❗✔️:     for (auto nbi2 = nbi3++; nbi2 != doneNbrs.end(); ++nbi2) {
        // RDKit❗✔️:       ang = computeAngle(d_eatoms[aid].loc, d_eatoms[*nbi1].loc,
        // RDKit❗✔️:                          d_eatoms[*nbi2].loc);
        // RDKit❗✔️:       auto nbrPair = std::make_pair((*nbi1), (*nbi2));
        // RDKit❗✔️:       anglePairs.emplace_back(ang, nbrPair);
        // RDKit❗✔️:     }
        // RDKit❗✔️:   }
        // RDKit❗✔️:   anglePairs.sort([](auto pr1, auto pr2) { return pr1.first < pr2.first; });
        // RDKit❗✔️:   auto winner = anglePairs.back();
        // RDKit❗✔️:   for (auto pr : boost::adaptors::reverse(anglePairs)) {
        // RDKit❗✔️:     if ((dp_mol->getRingInfo()->numAtomRings(pr.second.first) <= 1) &&
        // RDKit❗✔️:         (dp_mol->getRingInfo()->numAtomRings(pr.second.second) <= 1)) {
        // RDKit❗✔️:       winner = pr;
        // RDKit❗✔️:       break;
        // RDKit❗✔️:     }
        // RDKit❗✔️:   }
        // RDKit❗✔️:   auto winPair = winner.second;
        // RDKit❗✔️:   auto wnb1 = winPair.first;
        // RDKit❗✔️:   auto wnb2 = winPair.second;
        // RDKit❗✔️:   int nb2 = -1, nb1 = -1;
        // RDKit❗✔️:   for (auto anglePair : anglePairs) {
        // RDKit❗✔️:     auto nbrPair = anglePair.second;
        // RDKit❗✔️:     if (wnb1 == nbrPair.first) {
        // RDKit❗✔️:       nb2 = wnb1;
        // RDKit❗✔️:       nb1 = nbrPair.second;
        // RDKit❗✔️:       break;
        // RDKit❗✔️:     } else if (wnb1 == nbrPair.second) {
        // RDKit❗✔️:       nb2 = wnb1;
        // RDKit❗✔️:       nb1 = nbrPair.first;
        // RDKit❗✔️:       break;
        // RDKit❗✔️:     } else if (wnb2 == nbrPair.first) {
        // RDKit❗✔️:       nb2 = wnb2;
        // RDKit❗✔️:       nb1 = nbrPair.second;
        // RDKit❗✔️:       break;
        // RDKit❗✔️:     } else if (wnb2 == nbrPair.second) {
        // RDKit❗✔️:       nb2 = wnb2;
        // RDKit❗✔️:       nb1 = nbrPair.first;
        // RDKit❗✔️:       break;
        // RDKit❗✔️:     }
        // RDKit❗✔️:   }
        // RDKit❗✔️:   auto wAng = winner.first;
        // RDKit❗✔️:   d_eatoms[aid].rotDir = rotationDir(d_eatoms[aid].loc, d_eatoms[nb1].loc,
        // RDKit❗✔️:                                      d_eatoms[nb2].loc, wAng);
        // RDKit❗✔️:   d_eatoms[aid].nbr1 = nb1;
        // RDKit❗✔️:   d_eatoms[aid].nbr2 = nb2;
        // RDKit❗✔️:   d_eatoms[aid].angle = 2 * M_PI - wAng;
        // RDKit❗✔️: }
        if done_neighbors.len() < 3 {
            return Err(FragmentError::NotEnoughEmbeddedNeighbors {
                atom: aid,
                count: done_neighbors.len(),
            });
        }
        let center = self
            .atoms
            .get(&aid)
            .ok_or(FragmentError::AtomNotEmbedded { atom: aid })?
            .loc;
        let mut pairs = Vec::new();
        for (i, &first) in done_neighbors.iter().enumerate() {
            for &second in &done_neighbors[i..] {
                let a = self
                    .atoms
                    .get(&first)
                    .ok_or(FragmentError::AtomNotEmbedded { atom: first })?
                    .loc;
                let b = self
                    .atoms
                    .get(&second)
                    .ok_or(FragmentError::AtomNotEmbedded { atom: second })?
                    .loc;
                pairs.push((compute_angle(center, a, b)?, first, second));
            }
        }
        pairs.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(Ordering::Equal));
        let mut winner = *pairs.last().ok_or(FragmentError::InvalidAngle)?;
        for &pair in pairs.iter().rev() {
            if self.rings.num_atom_rings(AtomId::new(pair.1)) <= 1
                && self.rings.num_atom_rings(AtomId::new(pair.2)) <= 1
            {
                winner = pair;
                break;
            }
        }
        let (w_angle, win1, win2) = winner;
        let mut selected = None;
        for &(_, first, second) in &pairs {
            if win1 == first {
                selected = Some((second, win1));
            } else if win1 == second {
                selected = Some((first, win1));
            } else if win2 == first {
                selected = Some((second, win2));
            } else if win2 == second {
                selected = Some((first, win2));
            }
            if selected.is_some() {
                break;
            }
        }
        let (nb1, nb2) = selected.ok_or(FragmentError::InvalidAngle)?;
        let a = self.atoms[&nb1].loc;
        let b = self.atoms[&nb2].loc;
        let atom = self.atoms.get_mut(&aid).expect("checked above");
        atom.rot_dir = rotation_dir(center, a, b, w_angle);
        atom.nbr1 = Some(nb1);
        atom.nbr2 = Some(nb2);
        atom.angle = 2.0 * PI - w_angle;
        Ok(())
    }

    pub(crate) fn find_num_neigh(&self, point: Point2, radius: f64) -> usize {
        // RDKit❗✔️: int EmbeddedFrag::findNumNeigh(const RDGeom::Point2D &pt, double radius) {
        // RDKit❗✔️:   int res = 0;
        // RDKit❗✔️:   for (const auto &efi : d_eatoms) {
        // RDKit❗✔️:     const auto &rloc = efi.second.loc;
        // RDKit❗✔️:     if ((rloc - pt).length() < radius) {
        // RDKit❗✔️:       ++res;
        // RDKit❗✔️:     }
        // RDKit❗✔️:   }
        // RDKit❗✔️:   return res;
        // RDKit❗✔️: }
        self.atoms
            .values()
            .filter(|atom| {
                let delta = [atom.loc[0] - point[0], atom.loc[1] - point[1]];
                (delta[0] * delta[0] + delta[1] * delta[1]).sqrt() < radius
            })
            .count()
    }

    pub(crate) fn add_non_ring_atom(
        &mut self,
        aid: usize,
        to_aid: usize,
    ) -> Result<(), FragmentError> {
        // RDKit❗✔️: void EmbeddedFrag::addNonRingAtom(unsigned int aid, unsigned int toAid) {
        // RDKit❗✔️:   PRECONDITION(dp_mol, "");
        // RDKit❗✔️:   PRECONDITION(d_eatoms.find(aid) == d_eatoms.end(), "");
        // RDKit❗✔️:   PRECONDITION(d_eatoms.find(toAid) != d_eatoms.end(), "");
        // RDKit❗✔️:   if (d_eatoms[toAid].angle > 0.0) {
        // RDKit❗✔️:     addAtomToAtomWithAng(aid, toAid);
        // RDKit❗✔️:   } else {
        // RDKit❗✔️:     addAtomToAtomWithNoAng(aid, toAid);
        // RDKit❗✔️:   }
        // RDKit❗✔️:   d_eatoms[toAid].neighs.erase(std::remove(d_eatoms[toAid].neighs.begin(),
        // RDKit❗✔️:                                            d_eatoms[toAid].neighs.end(),
        // RDKit❗✔️:                                            static_cast<int>(aid)));
        // RDKit❗✔️:   this->updateNewNeighs(aid);
        // RDKit❗✔️: }
        if aid >= self.topology.atoms.len() || to_aid >= self.topology.atoms.len() {
            let atom = if aid >= self.topology.atoms.len() {
                aid
            } else {
                to_aid
            };
            return Err(FragmentError::AtomIndexOutOfRange {
                atom,
                atom_count: self.topology.atoms.len(),
            });
        }
        if self.atoms.contains_key(&aid) {
            return Err(FragmentError::AtomAlreadyEmbedded { atom: aid });
        }
        let angle = self
            .atoms
            .get(&to_aid)
            .ok_or(FragmentError::AtomNotEmbedded { atom: to_aid })?
            .angle;
        if angle > 0.0 {
            self.add_atom_to_atom_with_ang(aid, to_aid)?;
        } else {
            self.add_atom_to_atom_with_no_ang(aid, to_aid)?;
        }
        self.atoms
            .get_mut(&to_aid)
            .expect("checked above")
            .neighs
            .retain(|&neighbor| neighbor != aid);
        self.update_new_neighbors(aid)
    }

    fn add_atom_to_atom_with_ang(
        &mut self,
        aid: usize,
        to_aid: usize,
    ) -> Result<(), FragmentError> {
        // RDKit❗✔️: void EmbeddedFrag::addAtomToAtomWithAng(unsigned int aid, unsigned int toAid) {
        // RDKit❗✔️:   const auto &refAtom = d_eatoms[toAid];
        // RDKit❗✔️:   auto refLoc = refAtom.loc;
        // RDKit❗✔️:   RDGeom::Point2D origin(0.0, 0.0);
        // RDKit❗✔️:   PRECONDITION(refAtom.angle > 0.0, "");
        // RDKit❗✔️:   auto nnbr = refAtom.neighs.size();
        // RDKit❗✔️:   double remAngle = 2 * M_PI - refAtom.angle;
        // RDKit❗✔️:   auto currAngle = remAngle / (1 + nnbr);
        // RDKit❗✔️:   d_eatoms[toAid].angle += currAngle;
        // RDKit❗✔️:   const auto &nb1 = d_eatoms.at(refAtom.nbr1).loc;
        // RDKit❗✔️:   const auto &nb2 = d_eatoms.at(refAtom.nbr2).loc;
        // RDKit❗✔️:   if (d_eatoms[toAid].rotDir == 0) {
        // RDKit❗✔️:     d_eatoms[toAid].rotDir = rotationDir(refLoc, nb1, nb2, remAngle);
        // RDKit❗✔️:   }
        // RDKit❗✔️:   currAngle *= d_eatoms[toAid].rotDir;
        // RDKit❗✔️:   RDGeom::Transform2D rtrans;
        // RDKit❗✔️:   rtrans.SetTransform(refLoc, currAngle);
        // RDKit❗✔️:   auto currLoc = nb2;
        // RDKit❗✔️:   rtrans.TransformPoint(currLoc);
        // RDKit❗✔️:   if (fabs(remAngle) - M_PI < 1e-3) {
        // RDKit❗✔️:     auto currLoc2 = nb2;
        // RDKit❗✔️:     rtrans.SetTransform(refLoc, -currAngle);
        // RDKit❗✔️:     rtrans.TransformPoint(currLoc2);
        // RDKit❗✔️:     if (findNumNeigh(currLoc, 0.5) > findNumNeigh(currLoc2, 0.5)) {
        // RDKit❗✔️:       currLoc = currLoc2;
        // RDKit❗✔️:       currAngle *= -1;
        // RDKit❗✔️:     } else {
        // RDKit❗✔️:       rtrans.SetTransform(refLoc, currAngle);
        // RDKit❗✔️:     }
        // RDKit❗✔️:   }
        // RDKit❗✔️:   d_eatoms[toAid].nbr2 = aid;
        // RDKit❗✔️:   EmbeddedAtom eatm;
        // RDKit❗✔️:   eatm.aid = aid;
        // RDKit❗✔️:   eatm.loc = currLoc;
        // RDKit❗✔️:   eatm.nbr1 = toAid;
        // RDKit❗✔️:   eatm.angle = -1.0;
        // RDKit❗✔️:   auto tpt = currLoc - refLoc;
        // RDKit❗✔️:   RDGeom::Point2D norm(-tpt.y, tpt.x);
        // RDKit❗✔️:   auto tp1 = currLoc + norm;
        // RDKit❗✔️:   auto tp2 = currLoc - norm;
        // RDKit❗✔️:   auto nccw = findNumNeigh(tp1, NEIGH_RADIUS);
        // RDKit❗✔️:   auto ncw = findNumNeigh(tp2, NEIGH_RADIUS);
        // RDKit❗✔️:   norm.normalize();
        // RDKit❗✔️:   if (nccw < ncw) {
        // RDKit❗✔️:     eatm.normal = norm;
        // RDKit❗✔️:     eatm.ccw = false;
        // RDKit❗✔️:   } else {
        // RDKit❗✔️:     eatm.normal = (-norm);
        // RDKit❗✔️:     eatm.ccw = true;
        // RDKit❗✔️:   }
        // RDKit❗✔️:   d_eatoms[aid] = eatm;
        // RDKit❗✔️: }
        let reference = self
            .atoms
            .get(&to_aid)
            .ok_or(FragmentError::AtomNotEmbedded { atom: to_aid })?
            .clone();
        let first = reference
            .nbr1
            .ok_or(FragmentError::NotEnoughEmbeddedNeighbors {
                atom: to_aid,
                count: 0,
            })?;
        let second = reference
            .nbr2
            .ok_or(FragmentError::NotEnoughEmbeddedNeighbors {
                atom: to_aid,
                count: 1,
            })?;
        let nb1 = self
            .atoms
            .get(&first)
            .ok_or(FragmentError::AtomNotEmbedded { atom: first })?
            .loc;
        let nb2 = self
            .atoms
            .get(&second)
            .ok_or(FragmentError::AtomNotEmbedded { atom: second })?
            .loc;
        let rem_angle = 2.0 * PI - reference.angle;
        let mut current_angle = rem_angle / (1 + reference.neighs.len()) as f64;
        let mut rot_dir = reference.rot_dir;
        if rot_dir == 0 {
            rot_dir = rotation_dir(reference.loc, nb1, nb2, rem_angle);
        }
        current_angle *= f64::from(rot_dir);
        let mut current = Transform2D::around(reference.loc, current_angle).transform_point(nb2);
        if rem_angle.abs() - PI < 1e-3 {
            let other = Transform2D::around(reference.loc, -current_angle).transform_point(nb2);
            if self.find_num_neigh(current, 0.5) > self.find_num_neigh(other, 0.5) {
                current = other;
            }
        }
        let delta = [current[0] - reference.loc[0], current[1] - reference.loc[1]];
        let raw_normal = [-delta[1], delta[0]];
        let first_point = [current[0] + raw_normal[0], current[1] + raw_normal[1]];
        let second_point = [current[0] - raw_normal[0], current[1] - raw_normal[1]];
        let nccw = self.find_num_neigh(first_point, 2.5);
        let ncw = self.find_num_neigh(second_point, 2.5);
        let normal = normalize(raw_normal)?;
        let mut atom = EmbeddedAtom::at(aid, current);
        atom.nbr1 = Some(to_aid);
        if nccw < ncw {
            atom.normal = normal;
            atom.ccw = false;
        } else {
            atom.normal = [-normal[0], -normal[1]];
            atom.ccw = true;
        }
        let target = self.atoms.get_mut(&to_aid).expect("checked above");
        target.angle += rem_angle / (1 + reference.neighs.len()) as f64;
        target.rot_dir = rot_dir;
        target.nbr2 = Some(aid);
        self.atoms.insert(aid, atom);
        Ok(())
    }

    fn add_atom_to_atom_with_no_ang(
        &mut self,
        aid: usize,
        to_aid: usize,
    ) -> Result<(), FragmentError> {
        // RDKit❗✔️: void EmbeddedFrag::addAtomToAtomWithNoAng(unsigned int aid,
        // RDKit❗✔️:                                           unsigned int toAid) {
        // RDKit❗✔️:   PRECONDITION(dp_mol, "");
        // RDKit❗✔️:   const auto &refAtom = d_eatoms.at(toAid);
        // RDKit❗✔️:   PRECONDITION(refAtom.angle <= 0.0, "");
        // RDKit❗✔️:   const auto &refLoc = refAtom.loc;
        // RDKit❗✔️:   RDGeom::Point2D origin(0.0, 0.0);
        // RDKit❗✔️:   auto refAtomCCW = refAtom.ccw;
        // RDKit❗✔️:   auto currLoc = refAtom.normal;
        // RDKit❗✔️:   if (refAtom.CisTransNbr >= 0) {
        // RDKit❗✔️:     if (static_cast<unsigned int>(refAtom.CisTransNbr) != aid) {
        // RDKit❗✔️:       refAtomCCW = !refAtomCCW;
        // RDKit❗✔️:       currLoc *= -1.0;
        // RDKit❗✔️:     }
        // RDKit❗✔️:   }
        // RDKit❗✔️:   CHECK_INVARIANT(currLoc.lengthSq() > 1.0e-8, "");
        // RDKit❗✔️:   const auto atm = dp_mol->getAtomWithIdx(toAid);
        // RDKit❗✔️:   auto deg = getDepictDegree(atm);
        // RDKit❗✔️:   auto angle = computeSubAngle(deg, atm->getHybridization());
        // RDKit❗✔️:   bool flipNorm = false;
        // RDKit❗✔️:   if (d_eatoms[toAid].nbr1 >= 0) {
        // RDKit❗✔️:     d_eatoms[toAid].angle = angle;
        // RDKit❗✔️:     d_eatoms[toAid].nbr2 = aid;
        // RDKit❗✔️:   } else {
        // RDKit❗✔️:     auto norm = d_eatoms.at(toAid).normal;
        // RDKit❗✔️:     RDGeom::Transform2D rtrans;
        // RDKit❗✔️:     rtrans.SetTransform(origin, angle);
        // RDKit❗✔️:     rtrans.TransformPoint(norm);
        // RDKit❗✔️:     d_eatoms[toAid].normal = norm;
        // RDKit❗✔️:     d_eatoms[toAid].nbr1 = aid;
        // RDKit❗✔️:     flipNorm = true;
        // RDKit❗✔️:   }
        // RDKit❗✔️:   angle -= M_PI / 2;
        // RDKit❗✔️:   if (!refAtomCCW) {
        // RDKit❗✔️:     angle *= -1.0;
        // RDKit❗✔️:   }
        // RDKit❗✔️:   RDGeom::Transform2D trans;
        // RDKit❗✔️:   trans.SetTransform(origin, angle);
        // RDKit❗✔️:   trans.TransformPoint(currLoc);
        // RDKit❗✔️:   currLoc *= BOND_LEN;
        // RDKit❗✔️:   currLoc += refLoc;
        // RDKit❗✔️:   auto tpt = refLoc - currLoc;
        // RDKit❗✔️:   RDGeom::Point2D norm(-tpt.y, tpt.x);
        // RDKit❗✔️:   if (refAtomCCW ^ flipNorm) {
        // RDKit❗✔️:     norm *= -1.0;
        // RDKit❗✔️:   }
        // RDKit❗✔️:   norm.normalize();
        // RDKit❗✔️:   EmbeddedAtom eatm;
        // RDKit❗✔️:   eatm.loc = currLoc;
        // RDKit❗✔️:   eatm.normal = norm;
        // RDKit❗✔️:   eatm.nbr1 = toAid;
        // RDKit❗✔️:   eatm.angle = -1.0;
        // RDKit❗✔️:   eatm.ccw = (!refAtomCCW) ^ flipNorm;
        // RDKit❗✔️:   d_eatoms[aid] = eatm;
        // RDKit❗✔️: }
        let reference = self
            .atoms
            .get(&to_aid)
            .ok_or(FragmentError::AtomNotEmbedded { atom: to_aid })?
            .clone();
        let mut ccw = reference.ccw;
        let mut current = reference.normal;
        if reference
            .cis_trans_nbr
            .is_some_and(|neighbor| neighbor != aid)
        {
            ccw = !ccw;
            current = [-current[0], -current[1]];
        }
        if current[0] * current[0] + current[1] * current[1] <= 1e-8 {
            return Err(FragmentError::CoincidentPoints);
        }
        let degree = self.topology.adjacency.neighbors_of(to_aid).len();
        let mut angle = compute_sub_angle(degree, self.topology.atoms[to_aid].hybridization());
        let flip_normal = reference.nbr1.is_none();
        let new_ref_normal = if flip_normal {
            Some(Transform2D::around([0.0, 0.0], angle).transform_point(reference.normal))
        } else {
            None
        };
        angle -= PI / 2.0;
        if !ccw {
            angle *= -1.0;
        }
        current = Transform2D::around([0.0, 0.0], angle).transform_point(current);
        current = [
            current[0] * BOND_LEN + reference.loc[0],
            current[1] * BOND_LEN + reference.loc[1],
        ];
        let delta = [reference.loc[0] - current[0], reference.loc[1] - current[1]];
        let flip = ccw ^ flip_normal;
        let signed = if flip { -1.0 } else { 1.0 };
        let normal = normalize([-delta[1] * signed, delta[0] * signed])?;
        // The C++ default-constructed `eatm` does not assign its `aid` here.
        let mut atom = EmbeddedAtom::at(0, current);
        atom.normal = normal;
        atom.nbr1 = Some(to_aid);
        atom.ccw = (!ccw) ^ flip_normal;
        let target = self.atoms.get_mut(&to_aid).expect("checked above");
        if let Some(normal) = new_ref_normal {
            target.normal = normal;
            target.nbr1 = Some(aid);
        } else {
            target.angle = compute_sub_angle(degree, self.topology.atoms[to_aid].hybridization());
            target.nbr2 = Some(aid);
        }
        self.atoms.insert(aid, atom);
        Ok(())
    }
}

impl EmbeddedFrag<'_> {
    fn collect_flip_side(&self, end_aid: usize, begin_aid: usize) -> Vec<usize> {
        // RDKit❗✔️: void _recurseAtomOneSide(unsigned int endAid, unsigned int begAid,
        // RDKit❗✔️:                          const RDKit::ROMol *mol, RDKit::INT_VECT &flipAids) {
        // RDKit❗✔️:   PRECONDITION(mol, "");
        // RDKit❗✔️:   flipAids.push_back(endAid);
        // RDKit❗✔️:   for (auto nbr : mol->atomNeighbors(mol->getAtomWithIdx(endAid))) {
        // RDKit❗✔️:     if (nbr->getIdx() != begAid &&
        // RDKit❗✔️:         (std::find(flipAids.begin(), flipAids.end(),
        // RDKit❗✔️:                    static_cast<int>(nbr->getIdx())) == flipAids.end())) {
        // RDKit❗✔️:       _recurseAtomOneSide(nbr->getIdx(), begAid, mol, flipAids);
        // RDKit❗✔️:     }
        // RDKit❗✔️:   }
        // RDKit❗✔️:   return;
        // RDKit❗✔️: }
        // Behavior: source neighbor order and visited-vector test, with the
        // opposite bond endpoint excluded at every recursive level.
        // Complexity: source O(path traversal * visited-vector length).
        fn visit(topology: &TopologyBlock, end: usize, begin: usize, side: &mut Vec<usize>) {
            side.push(end);
            for neighbor in topology.adjacency.neighbors_of(end) {
                if neighbor.atom_index != begin && !side.contains(&neighbor.atom_index) {
                    visit(topology, neighbor.atom_index, begin, side);
                }
            }
        }
        let mut side = Vec::new();
        visit(self.topology, end_aid, begin_aid, &mut side);
        side
    }

    pub(crate) fn flip_about_bond(
        &mut self,
        bond_id: usize,
        flip_end: bool,
    ) -> Result<(), FragmentError> {
        // RDKit❗✔️: void EmbeddedFrag::flipAboutBond(unsigned int bondId, bool flipEnd) {
        // RDKit❗✔️:   PRECONDITION(dp_mol, "");
        // RDKit❗✔️:   PRECONDITION(bondId < dp_mol->getNumBonds(), "");
        // RDKit❗✔️:   // reflect all the atoms on one side of a bond using the bond as the mirror
        // RDKit❗✔️:   const auto bond = dp_mol->getBondWithIdx(bondId);
        // RDKit❗✔️:   // we should not be flip things around a ring bond
        // RDKit❗✔️:   CHECK_INVARIANT(!(dp_mol->getRingInfo()->numBondRings(bondId)), "");
        // RDKit❗✔️:   auto begAid = bond->getBeginAtomIdx();
        // RDKit❗✔️:   auto endAid = bond->getEndAtomIdx();
        // RDKit❗✔️:   if (!flipEnd) {
        // RDKit❗✔️:     std::swap(begAid, endAid);
        // RDKit❗✔️:   }
        // RDKit❗✔️:   const auto &begLoc = d_eatoms.at(begAid).loc;
        // RDKit❗✔️:   const auto &endLoc = d_eatoms.at(endAid).loc;
        // RDKit❗✔️:   // arbitrary choice here - find all atoms on one side of the bond
        // RDKit❗✔️:   // endAtom side - we will do this recursively
        // RDKit❗✔️:   RDKit::INT_VECT endSideAids;
        // RDKit❗✔️:   _recurseAtomOneSide(endAid, begAid, dp_mol, endSideAids);
        // RDKit❗✔️:   // look for fixed atoms in the fragment:
        // RDKit❗✔️:   unsigned int nAtomsFixed = 0;
        // RDKit❗✔️:   for (auto &d_eatom : d_eatoms) {
        // RDKit❗✔️:     if (d_eatom.second.df_fixed) {
        // RDKit❗✔️:       ++nAtomsFixed;
        // RDKit❗✔️:     }
        // RDKit❗✔️:   }
        // RDKit❗✔️:   // if there are fixed atoms, look at the atoms on the "end side"
        // RDKit❗✔️:   unsigned int nEndAtomsFixed = 0;
        // RDKit❗✔️:   if (nAtomsFixed) {
        // RDKit❗✔️:     for (auto endAtomId : endSideAids) {
        // RDKit❗✔️:       if (d_eatoms[endAtomId].df_fixed) {
        // RDKit❗✔️:         ++nEndAtomsFixed;
        // RDKit❗✔️:       }
        // RDKit❗✔️:     }
        // RDKit❗✔️:   }
        // RDKit❗✔️:   // now we have the molecule split into two groups of atoms
        // RDKit❗✔️:   // atom on the side of endAid and the rest.
        // RDKit❗✔️:   // we will flip the side that is smaller, assuming that there
        // RDKit❗✔️:   // are no fixed atoms there
        // RDKit❗✔️:   bool endSideFlip = true;
        // RDKit❗✔️:   if (nEndAtomsFixed) {
        // RDKit❗✔️:     endSideFlip = false;
        // RDKit❗✔️:     // there are fixed atoms on both sides, just return
        // RDKit❗✔️:     return;
        // RDKit❗✔️:   } else {
        // RDKit❗✔️:     auto nats = d_eatoms.size();
        // RDKit❗✔️:     auto nEndSide = endSideAids.size();
        // RDKit❗✔️:     if ((nats - nEndSide) < nEndSide) {
        // RDKit❗✔️:       endSideFlip = false;
        // RDKit❗✔️:     }
        // RDKit❗✔️:   }
        // RDKit❗✔️:   for (auto &d_eatom : d_eatoms) {
        // RDKit❗✔️:     const auto fii = std::find(endSideAids.begin(), endSideAids.end(),
        // RDKit❗✔️:                                static_cast<int>(d_eatom.first));
        // RDKit❗✔️:     if (endSideFlip ^ (fii == endSideAids.end())) {
        // RDKit❗✔️:       d_eatom.second.Reflect(begLoc, endLoc);
        // RDKit❗✔️:     }
        // RDKit❗✔️:   }
        // RDKit❗✔️: }
        // Behavior: source early return on any end-side fixed row and
        // smaller-side comparison; even a fixed row on the opposite side is
        // not separately protected by source if that side is smaller.
        // Complexity: one DFS, fixed-row scans and source-shaped visited tests.
        let bond = self
            .topology
            .bonds
            .get(bond_id)
            .ok_or(FragmentError::CollisionBondInvalid { bond: bond_id })?;
        if self.rings.num_bond_rings(bond.id()) != 0 {
            return Err(FragmentError::CollisionBondInvalid { bond: bond_id });
        }
        let (mut begin, mut end) = (bond.begin().index(), bond.end().index());
        if !flip_end {
            std::mem::swap(&mut begin, &mut end);
        }
        let begin_loc = self
            .atoms
            .get(&begin)
            .ok_or(FragmentError::AtomNotEmbedded { atom: begin })?
            .loc;
        let end_loc = self
            .atoms
            .get(&end)
            .ok_or(FragmentError::AtomNotEmbedded { atom: end })?
            .loc;
        let end_side = self.collect_flip_side(end, begin);
        if self.atoms.values().any(|atom| atom.fixed) {
            // RDKit❗✔️:   if (nAtomsFixed) {
            // RDKit❗✔️:     for (auto endAtomId : endSideAids) {
            // RDKit❗✔️:       if (d_eatoms[endAtomId].df_fixed) {
            // RDKit❗✔️:         ++nEndAtomsFixed;
            // RDKit❗✔️:       }
            // RDKit❗✔️:     }
            // RDKit❗✔️:   }
            // Behavior: operator[] inserts source-default rows for missing
            // recursive-side atoms before the fixed-side decision.
            // Complexity: one ordered-map lookup/insertion per recursive row,
            // matching std::map::operator[].
            let mut end_has_fixed = false;
            for &atom_id in &end_side {
                let atom = self
                    .atoms
                    .entry(atom_id)
                    .or_insert_with(EmbeddedAtom::source_default);
                end_has_fixed |= atom.fixed;
            }
            if end_has_fixed {
                return Ok(());
            }
        }
        // `size_t` subtraction is unsigned in the source. This is observable
        // for a partial fragment with no fixed rows and a larger recursive
        // whole-molecule side; wrapping avoids a debug-only invented panic.
        let end_side_flip = self.atoms.len().wrapping_sub(end_side.len()) >= end_side.len();
        for (&aid, atom) in &mut self.atoms {
            if end_side_flip ^ !end_side.contains(&aid) {
                atom.reflect(begin_loc, end_loc);
            }
        }
        Ok(())
    }

    fn all_rotatable_bonds(&self) -> Vec<usize> {
        // RDKit❗✔️: RDKit::INT_VECT getAllRotatableBonds(const RDKit::ROMol &mol) {
        // RDKit❗✔️:   RDKit::INT_VECT res;
        // RDKit❗✔️:   for (const auto bond : mol.bonds()) {
        // RDKit❗✔️:     int bid = bond->getIdx();
        // RDKit❗✔️:     if ((bond->getStereo() <= RDKit::Bond::STEREOANY) &&
        // RDKit❗✔️:         (!(mol.getRingInfo()->numBondRings(bid)))) {
        // RDKit❗✔️:       res.push_back(bid);
        // RDKit❗✔️:     }
        // RDKit❗✔️:   }
        // RDKit❗✔️:   return res;
        // RDKit❗✔️: }
        // Behavior: the candidate domain is the whole molecule in bond-table
        // order. It is deliberately not filtered to the current fragment.
        // Complexity: one topology-bond scan and one output allocation, as in
        // the source vector construction.
        self.topology
            .bonds
            .iter()
            .filter(|bond| {
                bond.stereo().rdkit_code() <= BondStereo::Any.rdkit_code()
                    && self.rings.num_bond_rings(bond.id()) == 0
            })
            .map(|bond| bond.id().index())
            .collect()
    }

    fn degree_four_candidates(&self) -> Vec<DegreeFourCandidate> {
        // RDKit❗✔️:   if (permuteDeg4Nodes) {
        // RDKit❗✔️:     for (const auto atom : dp_mol->atoms()) {
        // RDKit❗✔️:       auto caid = atom->getIdx();
        // RDKit❗✔️:       if ((getDepictDegree(atom) == 4) &&
        // RDKit❗✔️:           (!(dp_mol->getRingInfo()->numAtomRings(caid)))) {
        // RDKit❗✔️:         RDKit::INT_VECT aids, bids;
        // RDKit❗✔️:         getNbrAtomAndBondIds(caid, dp_mol, aids, bids);
        // RDKit❗✔️:         bool allin = true;
        // RDKit❗✔️:         for (auto aid : aids) {
        // RDKit❗✔️:           auto nbrIter = d_eatoms.find(aid);
        // RDKit❗✔️:           if (nbrIter == d_eatoms.end() || nbrIter->second.df_fixed) {
        // RDKit❗✔️:             allin = false;
        // RDKit❗✔️:             break;
        // RDKit❗✔️:           }
        // RDKit❗✔️:         }
        // RDKit❗✔️:         if (allin) {
        // RDKit❗✔️:           deg4nodes.push_back(caid);
        // RDKit❗✔️:           deg4NbrBids.push_back(bids);
        // RDKit❗✔️:           deg4NbrAids.push_back(aids);
        // RDKit❗✔️:         }
        // RDKit❗✔️:       }
        // RDKit❗✔️:     }
        // RDKit❗✔️:   }
        // RDKit❗✔️: void getNbrAtomAndBondIds(unsigned int aid, const RDKit::ROMol *mol,
        // RDKit❗✔️:                           RDKit::INT_VECT &aids, RDKit::INT_VECT &bids) {
        // RDKit❗✔️:   for (auto nbr : mol->atomNeighbors(mol->getAtomWithIdx(aid))) {
        // RDKit❗✔️:     auto bi = mol->getBondBetweenAtoms(aid, nbr->getIdx())->getIdx();
        // RDKit❗✔️:     aids.push_back(nbr->getIdx());
        // RDKit❗✔️:     bids.push_back(bi);
        // RDKit❗✔️:   }
        // RDKit❗✔️: }
        // Behavior: atom and adjacency order are canonical table order. Only
        // neighbor rows, not the centre row, participate in the source all-in
        // and fixed test.
        // Complexity: one atom scan and one four-row adjacency scan per
        // eligible centre, matching the source shape.
        let mut candidates = Vec::new();
        for centre in 0..self.topology.atoms.len() {
            let neighbors = self.topology.adjacency.neighbors_of(centre);
            if neighbors.len() != 4 || self.rings.num_atom_rings(AtomId::new(centre)) != 0 {
                continue;
            }
            if neighbors.iter().any(|neighbor| {
                self.atoms
                    .get(&neighbor.atom_index)
                    .is_none_or(|atom| atom.fixed)
            }) {
                continue;
            }
            candidates.push(DegreeFourCandidate {
                centre,
                neighbor_atoms: neighbors.iter().map(|row| row.atom_index).collect(),
                neighbor_bonds: neighbors.iter().map(|row| row.bond.index()).collect(),
            });
        }
        candidates
    }

    fn degree_four_bond_pairs(
        &self,
        candidate: &DegreeFourCandidate,
    ) -> Result<[(usize, usize); 2], FragmentError> {
        // RDKit❗✔️:   double dp1 = nbrPts[0].dotProduct(nbrPts[1]);
        // RDKit❗✔️:   if (fabs(dp1) < 1.e-3) {
        // RDKit❗✔️:     INT_PAIR p1(nbrBids[0], nbrBids[1]);
        // RDKit❗✔️:     res.push_back(p1);
        // RDKit❗✔️:     double dp2 = nbrPts[0].dotProduct(nbrPts[2]);
        // RDKit❗✔️:     if (fabs(dp2) < 1.e-3) {
        // RDKit❗✔️:       INT_PAIR p2(nbrBids[0], nbrBids[2]);
        // RDKit❗✔️:       res.push_back(p2);
        // RDKit❗✔️:     } else {
        // RDKit❗✔️:       INT_PAIR p2(nbrBids[0], nbrBids[3]);
        // RDKit❗✔️:       res.push_back(p2);
        // RDKit❗✔️:     }
        // RDKit❗✔️:     return res;
        // RDKit❗✔️:   } else {
        // RDKit❗✔️:     INT_PAIR p1(nbrBids[0], nbrBids[2]);
        // RDKit❗✔️:     res.push_back(p1);
        // RDKit❗✔️:     INT_PAIR p2(nbrBids[0], nbrBids[3]);
        // RDKit❗✔️:     res.push_back(p2);
        // RDKit❗✔️:     return res;
        // RDKit❗✔️:   }
        // Behavior: exactly the first-neighbor dot-product decision tree and
        // strict 1e-3 threshold from the source helper.
        // Complexity: four map lookups and two constant-time dot products.
        let centre = self
            .atoms
            .get(&candidate.centre)
            .ok_or(FragmentError::AtomNotEmbedded {
                atom: candidate.centre,
            })?
            .loc;
        let mut vectors = [[0.0; 2]; 4];
        for (slot, &atom_id) in candidate.neighbor_atoms.iter().enumerate() {
            let loc = self
                .atoms
                .get(&atom_id)
                .ok_or(FragmentError::AtomNotEmbedded { atom: atom_id })?
                .loc;
            vectors[slot] = [loc[0] - centre[0], loc[1] - centre[1]];
        }
        let dot01 = vectors[0][0] * vectors[1][0] + vectors[0][1] * vectors[1][1];
        let pairs = if dot01.abs() < 1.0e-3 {
            let dot02 = vectors[0][0] * vectors[2][0] + vectors[0][1] * vectors[2][1];
            [
                (candidate.neighbor_bonds[0], candidate.neighbor_bonds[1]),
                (
                    candidate.neighbor_bonds[0],
                    if dot02.abs() < 1.0e-3 {
                        candidate.neighbor_bonds[2]
                    } else {
                        candidate.neighbor_bonds[3]
                    },
                ),
            ]
        } else {
            [
                (candidate.neighbor_bonds[0], candidate.neighbor_bonds[2]),
                (candidate.neighbor_bonds[0], candidate.neighbor_bonds[3]),
            ]
        };
        Ok(pairs)
    }

    pub(crate) fn permute_bonds(
        &mut self,
        centre: usize,
        first: usize,
        second: usize,
    ) -> Result<(), FragmentError> {
        // RDKit❗✔️: void EmbeddedFrag::permuteBonds(unsigned int aid, unsigned int aid1,
        // RDKit❗✔️:                                 unsigned int aid2) {
        // RDKit❗✔️:   auto rl1 = d_eatoms.at(aid).loc;
        // RDKit❗✔️:   auto rl2 = d_eatoms.at(aid1).loc + d_eatoms.at(aid2).loc;
        // RDKit❗✔️:   rl2 *= 0.5;
        // RDKit❗✔️:   RDKit::INT_VECT fragA, fragB;
        // RDKit❗✔️:   _recurseAtomOneSide(aid1, aid, dp_mol, fragA);
        // RDKit❗✔️:   _recurseAtomOneSide(aid2, aid, dp_mol, fragB);
        // RDKit❗✔️:   for (auto fi : fragA) {
        // RDKit❗✔️:     d_eatoms[fi].Reflect(rl1, rl2);
        // RDKit❗✔️:   }
        // RDKit❗✔️:   for (auto fi : fragB) {
        // RDKit❗✔️:     d_eatoms[fi].Reflect(rl1, rl2);
        // RDKit❗✔️:   }
        // RDKit❗✔️: }
        // Behavior: both recursive sides are reflected in source traversal
        // order around the centre-to-neighbor-midpoint bisector.
        // Complexity: two source-shaped DFS traversals and one reflection per
        // visited row; no topology clone.
        let centre_loc = self
            .atoms
            .get(&centre)
            .ok_or(FragmentError::AtomNotEmbedded { atom: centre })?
            .loc;
        let first_loc = self
            .atoms
            .get(&first)
            .ok_or(FragmentError::AtomNotEmbedded { atom: first })?
            .loc;
        let second_loc = self
            .atoms
            .get(&second)
            .ok_or(FragmentError::AtomNotEmbedded { atom: second })?
            .loc;
        let midpoint = [
            (first_loc[0] + second_loc[0]) * 0.5,
            (first_loc[1] + second_loc[1]) * 0.5,
        ];
        let first_side = self.collect_flip_side(first, centre);
        let second_side = self.collect_flip_side(second, centre);
        // RDKit❗✔️:   for (auto fi : fragA) {
        // RDKit❗✔️:     d_eatoms[fi].Reflect(rl1, rl2);
        // RDKit❗✔️:   }
        // RDKit❗✔️:   for (auto fi : fragB) {
        // RDKit❗✔️:     d_eatoms[fi].Reflect(rl1, rl2);
        // RDKit❗✔️:   }
        // Behavior: source operator[] inserts a default EmbeddedAtom for a
        // recursive-side row absent from a partial fragment before reflecting
        // it. Initial centre/neighbor lookups above remain source `.at()`
        // failure sites and are deliberately not converted to insertion.
        // Complexity: one ordered-map lookup/insertion and reflection per DFS
        // row, matching the two source loops.
        for atom_id in first_side.into_iter().chain(second_side) {
            self.atoms
                .entry(atom_id)
                .or_insert_with(EmbeddedAtom::source_default)
                .reflect(centre_loc, midpoint);
        }
        Ok(())
    }

    fn sampling_cost(
        &self,
        target_distances: Option<&[f64]>,
        mimic_distance_weight: f64,
    ) -> Result<f64, FragmentError> {
        // RDKit❗✔️: auto na = dp_mol->getNumAtoms();
        // RDKit❗✔️: if (na < 2) {
        // RDKit❗✔️:   return 0;
        // RDKit❗✔️: }
        // RDKit❗✔️: auto dsize = na * (na - 1) / 2;
        // RDKit❗✔️: auto *ddata2D = new double[dsize];
        // RDKit❗✔️: DOUBLE_SMART_PTR dmat2D(ddata2D);
        // RDKit❗✔️: this->computeDistMat(dmat2D);
        // RDKit❗✔️: double res1 = 0.0;
        // RDKit❗✔️: double res2 = 0.0;
        // RDKit❗✔️: for (auto i = 0u; i < dsize; ++i) {
        // RDKit❗✔️:   auto d = ddata2D[i];
        // RDKit❗✔️:   auto d2 = d * d;
        // RDKit❗✔️:   if (d2 > 1.e-3) {
        // RDKit❗✔️:     res1 += 1.0 / d2;
        // RDKit❗✔️:   } else {
        // RDKit❗✔️:     res1 += 1000.0;
        // RDKit❗✔️:   }
        // RDKit❗✔️:   if (ddata && (ddata[i] >= 0.0)) {
        // RDKit❗✔️:     auto dd = d - ddata[i];
        // RDKit❗✔️:     res2 += dd * dd;
        // RDKit❗✔️:   }
        // RDKit❗✔️: }
        // RDKit❗✔️: auto wt = mimicDmatWt;
        // RDKit❗✔️: if (wt > 1.0) {
        // RDKit❗✔️:   wt = 1.0;
        // RDKit❗✔️: } else if (wt < 0.0) {
        // RDKit❗✔️:   wt = 0.0;
        // RDKit❗✔️: }
        // RDKit❗✔️: return ((1.0 - wt) * res1) + (wt * res2);
        // RDKit❗✔️: void EmbeddedFrag::computeDistMat(DOUBLE_SMART_PTR &dmat) {
        // RDKit❗✔️:   for (auto efi = d_eatoms.begin(); efi != d_eatoms.end(); ++efi) {
        // RDKit❗✔️:     auto pti = efi->second.loc;
        // RDKit❗✔️:     auto ai = efi->first;
        // RDKit❗✔️:     for (auto efj = d_eatoms.begin(); efj != efi; ++efj) {
        // RDKit❗✔️:       auto ptj = efj->second.loc;
        // RDKit❗✔️:       auto aj = efj->first;
        // RDKit❗✔️:       ptj -= pti;
        // RDKit❗✔️:       if (ai < aj) {
        // RDKit❗✔️:         std::swap(ai, aj);
        // RDKit❗✔️:       }
        // RDKit❗✔️:       dmatPtr[(ai * (ai - 1) / 2) + aj] = ptj.length();
        // RDKit❗✔️:     }
        // RDKit❗✔️:   }
        // RDKit❗✔️: }
        // Behavior: complete fragments reproduce the triangular fill/read and
        // accumulation order. Option explicitly tracks unwritten partial rows
        // rather than reading C++ uninitialized storage or inventing a value.
        // Complexity: O(E^2 + V^2) time and O(V^2) storage, matching source.
        let atom_count = self.topology.atoms.len();
        if atom_count < 2 {
            return Ok(0.0);
        }
        let distance_count = atom_count * (atom_count - 1) / 2;
        if target_distances.is_some_and(|values| values.len() < distance_count) {
            return Err(FragmentError::MismatchedTopology);
        }
        let mut distances = vec![None; distance_count];
        let embedded: Vec<_> = self.atoms.keys().copied().collect();
        for (outer, &first_id) in embedded.iter().enumerate() {
            let first = &self.atoms[&first_id];
            for &second_id in &embedded[..outer] {
                let second = &self.atoms[&second_id];
                let (high, low) = if first_id > second_id {
                    (first_id, second_id)
                } else {
                    (second_id, first_id)
                };
                let dx = second.loc[0] - first.loc[0];
                let dy = second.loc[1] - first.loc[1];
                distances[high * (high - 1) / 2 + low] = Some(dx.hypot(dy));
            }
        }
        let mut density = 0.0;
        let mut mimic = 0.0;
        for (index, value) in distances.into_iter().enumerate() {
            let distance = value.ok_or_else(|| {
                let mut high = 1;
                while high * (high - 1) / 2 + high <= index {
                    high += 1;
                }
                let low = index - high * (high - 1) / 2;
                FragmentError::UndefinedSamplingDistance {
                    first: high,
                    second: low,
                }
            })?;
            let squared = distance * distance;
            density += if squared > 1.0e-3 {
                1.0 / squared
            } else {
                1000.0
            };
            if let Some(target) = target_distances
                .and_then(|values| values.get(index))
                .copied()
                .filter(|value| *value >= 0.0)
            {
                let delta = distance - target;
                mimic += delta * delta;
            }
        }
        let weight = mimic_distance_weight.clamp(0.0, 1.0);
        Ok((1.0 - weight) * density + weight * mimic)
    }

    pub(crate) fn random_sample_flips_and_permutations(
        &mut self,
        bonds_per_sample: u32,
        samples: u32,
        seed: i32,
        target_distances: Option<&[f64]>,
        mimic_distance_weight: f64,
        permute_degree_four: bool,
    ) -> Result<(), FragmentError> {
        // RDKit❗✔️: const auto &rotBonds = getAllRotatableBonds(*dp_mol);
        // RDKit❗✔️: auto nb = rotBonds.size();
        // RDKit❗✔️: unsigned int nt = nb + nd4;
        // RDKit❗✔️: unsigned int nPerSample = std::min(nt, nBondsPerSample);
        // RDKit❗✔️: auto &generator = RDKit::getRandomGenerator();
        // RDKit❗✔️: if (seed > 0) {
        // RDKit❗✔️:   generator.seed(seed);
        // RDKit❗✔️: }
        // RDKit❗✔️: RDKit::uniform_int dist(0, nt - 1);
        // RDKit❗✔️: RDKit::int_source_type intRandomSrc(generator, dist);
        // RDKit❗✔️: RDGeom::INT_POINT2D_MAP bestCrdMap;
        // RDKit❗✔️: auto bestDens = this->mimicDistMatAndDensityCostFunc(dmat, mimicDmatWt);
        // RDKit❗✔️: for (const auto &efi : d_eatoms) {
        // RDKit❗✔️:   bestCrdMap[efi.first] = efi.second.loc;
        // RDKit❗✔️: }
        // RDKit❗✔️: for (auto si = 0u; si < nSamples; ++si) {
        // RDKit❗✔️:   for (auto fi = 0u; fi < nPerSample; ++fi) {
        // RDKit❗✔️:     unsigned int ri = intRandomSrc();
        // RDKit❗✔️:     if (ri < nb) {
        // RDKit❗✔️:       this->flipAboutBond(rotBonds.at(ri));
        // RDKit❗✔️:     } else {
        // RDKit❗✔️:       unsigned int d4i = ri - nb;
        // RDKit❗✔️:       auto ai = deg4nodes.at(d4i);
        // RDKit❗✔️:       auto bndPairs = findBondsPairsToPermuteDeg4(
        // RDKit❗✔️:           d_eatoms.at(ai).loc, deg4NbrBids.at(d4i), nbrLocs);
        // RDKit❗✔️:       auto rval = RDKit::getRandomVal();
        // RDKit❗✔️:       unsigned int fbi = 0;
        // RDKit❗✔️:       if (rval > 0.5) {
        // RDKit❗✔️:         fbi = 1;
        // RDKit❗✔️:       }
        // RDKit❗✔️:       auto aid1 = dp_mol->getBondWithIdx(bndPairs.at(fbi).first)->getOtherAtomIdx(ai);
        // RDKit❗✔️:       auto aid2 = dp_mol->getBondWithIdx(bndPairs.at(fbi).second)->getOtherAtomIdx(ai);
        // RDKit❗✔️:       this->permuteBonds(ai, aid1, aid2);
        // RDKit❗✔️:     }
        // RDKit❗✔️:   }
        // RDKit❗✔️:   auto density = this->mimicDistMatAndDensityCostFunc(dmat, mimicDmatWt);
        // RDKit❗✔️:   if (bestDens - density > 1e-4) {
        // RDKit❗✔️:     bestDens = density;
        // RDKit❗✔️:     for (const auto &efi : d_eatoms) {
        // RDKit❗✔️:       bestCrdMap[efi.first] = efi.second.loc;
        // RDKit❗✔️:     }
        // RDKit❗✔️:   }
        // RDKit❗✔️: }
        // RDKit❗✔️: for (auto &efi : d_eatoms) {
        // RDKit❗✔️:   efi.second.loc = bestCrdMap.at(efi.first);
        // RDKit❗✔️: }
        // Behavior: whole-molecule candidate order, one shared Boost stream,
        // cumulative per-sample moves, strict improvement threshold and final
        // coordinate-only restore follow source. Partial cost remains an
        // explicit undefined-source boundary in sampling_cost.
        // Complexity: source-equivalent candidate scans, trial transforms and
        // quadratic cost computations; best coordinates allocate one map.
        let rotatable_bonds = self.all_rotatable_bonds();
        let degree_four = if permute_degree_four {
            self.degree_four_candidates()
        } else {
            Vec::new()
        };
        let rotatable_count = rotatable_bonds.len();
        let candidate_count = rotatable_count + degree_four.len();
        let per_sample = candidate_count.min(bonds_per_sample as usize);
        let mut generator = DEPICT_RANDOM.lock().expect("depiction RNG mutex poisoned");
        generator.reset_if_positive(seed);
        let mut best_cost = self.sampling_cost(target_distances, mimic_distance_weight)?;
        let mut best_coordinates: PointMap = self
            .atoms
            .iter()
            .map(|(&atom_id, atom)| (atom_id, atom.loc))
            .collect();
        for _ in 0..samples {
            for _ in 0..per_sample {
                let selected = generator.next_index(candidate_count);
                if selected < rotatable_count {
                    self.flip_about_bond(rotatable_bonds[selected], true)?;
                } else {
                    let candidate = &degree_four[selected - rotatable_count];
                    let pairs = self.degree_four_bond_pairs(candidate)?;
                    let pair_index = usize::from(generator.next_real() > 0.5);
                    let (first_bond, second_bond) = pairs[pair_index];
                    let other_endpoint = |bond_id: usize| {
                        let bond = &self.topology.bonds[bond_id];
                        if bond.begin().index() == candidate.centre {
                            Ok(bond.end().index())
                        } else if bond.end().index() == candidate.centre {
                            Ok(bond.begin().index())
                        } else {
                            Err(FragmentError::CollisionBondInvalid { bond: bond_id })
                        }
                    };
                    let first = other_endpoint(first_bond)?;
                    let second = other_endpoint(second_bond)?;
                    self.permute_bonds(candidate.centre, first, second)?;
                }
            }
            let cost = self.sampling_cost(target_distances, mimic_distance_weight)?;
            if best_cost - cost > 1.0e-4 {
                best_cost = cost;
                for (&atom_id, atom) in &self.atoms {
                    best_coordinates.insert(atom_id, atom.loc);
                }
            }
        }
        for (&atom_id, atom) in &mut self.atoms {
            atom.loc = best_coordinates[&atom_id];
        }
        Ok(())
    }

    pub(crate) fn remove_collisions_bond_flip(&mut self) -> Result<(), FragmentError> {
        // RDKit❗✔️: void EmbeddedFrag::removeCollisionsBondFlip() {
        // RDKit❗✔️:   // try to remove collisions in a structure by flipping rotatable bonds along
        // RDKit❗✔️:   // the shortest path between the colliding atoms. we will limit the number of
        // RDKit❗✔️:   // times we are going to do this since we may fall into spiral where removing
        // RDKit❗✔️:   // a collision may create a new one
        // RDKit❗✔️:   auto dmat = RDKit::MolOps::getDistanceMat(*dp_mol);
        // RDKit❗✔️:   auto colls = this->findCollisions(dmat);
        // RDKit❗✔️:   std::map<int, unsigned int> doneBonds;
        // RDKit❗✔️:   unsigned int iter = 0;
        // RDKit❗✔️:   while (iter < MAX_COLL_ITERS && colls.size()) {
        // RDKit❗✔️:     auto ncols = colls.size();
        // RDKit❗✔️:     if (ncols > 0) {
        // RDKit❗✔️:       // we have a collision
        // RDKit❗✔️:       auto cAids = colls[0];
        // RDKit❗✔️:       auto rotBonds = getRotatableBonds(*dp_mol, cAids.first, cAids.second);
        // RDKit❗✔️:       auto prevDensity = this->totalDensity();
        // RDKit❗✔️:       for (auto ri : rotBonds) {
        // RDKit❗✔️:         auto doneBondsRiIt = doneBonds.find(ri);
        // RDKit❗✔️:         if ((doneBondsRiIt == doneBonds.end()) ||
        // RDKit❗✔️:             (doneBondsRiIt->second < NUM_BONDS_FLIPS)) {
        // RDKit❗✔️:           if (doneBondsRiIt == doneBonds.end()) {
        // RDKit❗✔️:             doneBonds[ri] = 1;
        // RDKit❗✔️:           } else {
        // RDKit❗✔️:             doneBondsRiIt->second += 1;
        // RDKit❗✔️:           }
        // RDKit❗✔️:           flipAboutBond(ri);
        // RDKit❗✔️:           colls = this->findCollisions(dmat);
        // RDKit❗✔️:           auto newDensity = this->totalDensity();
        // RDKit❗✔️:           if (colls.size() < ncols) {
        // RDKit❗✔️:             doneBonds[ri] = NUM_BONDS_FLIPS;  // lock this rotatable bond
        // RDKit❗✔️:             break;
        // RDKit❗✔️:           } else if (colls.size() == ncols && newDensity < prevDensity) {
        // RDKit❗✔️:             break;
        // RDKit❗✔️:           } else {
        // RDKit❗✔️:             // we made the wrong move earlier - reject the flip move it back
        // RDKit❗✔️:             flipAboutBond(ri);
        // RDKit❗✔️:             colls = this->findCollisions(dmat);
        // RDKit❗✔️:             // and try the other end:
        // RDKit❗✔️:             flipAboutBond(ri, false);
        // RDKit❗✔️:             colls = this->findCollisions(dmat);
        // RDKit❗✔️:             newDensity = this->totalDensity();
        // RDKit❗✔️:             if (colls.size() < ncols) {
        // RDKit❗✔️:               doneBonds[ri] = NUM_BONDS_FLIPS;  // lock this rotatable bond
        // RDKit❗✔️:               break;
        // RDKit❗✔️:             } else if (colls.size() == ncols && newDensity < prevDensity) {
        // RDKit❗✔️:               break;
        // RDKit❗✔️:             } else {
        // RDKit❗✔️:               flipAboutBond(ri, false);
        // RDKit❗✔️:               colls = this->findCollisions(dmat);
        // RDKit❗✔️:             }
        // RDKit❗✔️:           }
        // RDKit❗✔️:         }
        // RDKit❗✔️:       }
        // RDKit❗✔️:     }
        // RDKit❗✔️:     ++iter;
        // RDKit❗✔️:   }
        // RDKit❗✔️: }
        // Behavior: source first-collision selection, three flips per bond,
        // 15 global iterations, two trial orientations and exact rejection.
        // Complexity: each trial recomputes collision/density as source.
        let distance = self.collision_distance_matrix()?;
        let mut collisions = self.find_collisions(&distance, true);
        let mut done_bonds = BTreeMap::<usize, usize>::new();
        let mut iteration = 0;
        while iteration < 15 && !collisions.is_empty() {
            let old_count = collisions.len();
            let (first, second) = collisions[0];
            let rotatable = self.rotatable_bonds_on_shortest_path(first, second)?;
            let old_density = self.total_density();
            for bond in rotatable {
                if done_bonds.get(&bond).is_some_and(|&count| count >= 3) {
                    continue;
                }
                *done_bonds.entry(bond).or_default() += 1;
                self.flip_about_bond(bond, true)?;
                collisions = self.find_collisions(&distance, true);
                let new_density = self.total_density();
                if collisions.len() < old_count {
                    done_bonds.insert(bond, 3);
                    break;
                }
                if collisions.len() == old_count && new_density < old_density {
                    break;
                }
                self.flip_about_bond(bond, true)?;
                collisions = self.find_collisions(&distance, true);
                self.flip_about_bond(bond, false)?;
                collisions = self.find_collisions(&distance, true);
                let new_density = self.total_density();
                if collisions.len() < old_count {
                    done_bonds.insert(bond, 3);
                    break;
                }
                if collisions.len() == old_count && new_density < old_density {
                    break;
                }
                self.flip_about_bond(bond, false)?;
                collisions = self.find_collisions(&distance, true);
            }
            iteration += 1;
        }
        Ok(())
    }

    pub(crate) fn collision_distance_matrix(&self) -> Result<DenseMatrix, FragmentError> {
        // RDKit❗✔️:   auto dmat = RDKit::MolOps::getDistanceMat(*dp_mol);
        // Behavior: the default full topological distance matrix is the
        // source graph-distance input, independent of current XY positions.
        // Complexity: reuse the canonical core Floyd-Warshall owner.
        Ok(topological_distance_matrix(
            self.topology,
            &TopologicalDistanceMatrixParams::default(),
        )?)
    }

    fn degree_one_neighbor(&self, aid: usize) -> Result<usize, FragmentError> {
        // RDKit❗✔️: unsigned int _findDeg1Neighbor(const RDKit::ROMol *mol, unsigned int aid) {
        // RDKit❗✔️:   PRECONDITION(mol, "");
        // RDKit❗✔️:   auto deg = getDepictDegree(mol->getAtomWithIdx(aid));
        // RDKit❗✔️:   CHECK_INVARIANT(deg == 1, "");
        // RDKit❗✔️:   return *mol->getAtomNeighbors(mol->getAtomWithIdx(aid)).first;
        // RDKit❗✔️: }
        // Behavior: getDepictDegree is atom graph degree in this pinned source.
        // Complexity: direct first adjacency row, no graph traversal.
        let neighbors = self.topology.adjacency.neighbors_of(aid);
        if neighbors.len() != 1 {
            return Err(FragmentError::NotEnoughEmbeddedNeighbors {
                atom: aid,
                count: neighbors.len(),
            });
        }
        Ok(neighbors[0].atom_index)
    }

    fn closest_neighbor(&self, distance: &DenseMatrix, aid1: usize, aid2: usize) -> usize {
        // RDKit❗✔️: unsigned int _findClosestNeighbor(const RDKit::ROMol *mol, const double *dmat,
        // RDKit❗✔️:                                   unsigned int aid1, unsigned int aid2) {
        // RDKit❗✔️:   PRECONDITION(mol, "");
        // RDKit❗✔️:   unsigned int res = 0;
        // RDKit❗✔️:   double mdist = 1.e8;
        // RDKit❗✔️:   auto naid = aid1 * (mol->getNumAtoms());
        // RDKit❗✔️:   for (const auto nbr : mol->atomNeighbors(mol->getAtomWithIdx(aid2))) {
        // RDKit❗✔️:     auto d = dmat[naid + nbr->getIdx()];
        // RDKit❗✔️:     if (d < mdist) {
        // RDKit❗✔️:       mdist = d;
        // RDKit❗✔️:       res = nbr->getIdx();
        // RDKit❗✔️:     }
        // RDKit❗✔️:   }
        // RDKit❗✔️:   return res;
        // RDKit❗✔️: }
        // Behavior: retain source default zero and strict first-minimum tie.
        // Complexity: one adjacency scan and O(1) matrix access per neighbor.
        let mut result = 0;
        let mut minimum = 1.0e8;
        for neighbor in self.topology.adjacency.neighbors_of(aid2) {
            let value = distance
                .get(aid1, neighbor.atom_index)
                .expect("full matrix");
            if value < minimum {
                minimum = value;
                result = neighbor.atom_index;
            }
        }
        result
    }

    pub(crate) fn open_angles(
        &mut self,
        distance: &DenseMatrix,
        aid1: usize,
        aid2: usize,
    ) -> Result<(), FragmentError> {
        // RDKit❗✔️: void EmbeddedFrag::openAngles(const double *dmat, unsigned int aid1,
        // RDKit❗✔️:                               unsigned int aid2) {
        // RDKit❗✔️:   // Assuming that either aid1, and/or aid2 are degree 1 atoms, we will open up
        // RDKit❗✔️:   // the angles
        // RDKit❗✔️:   //
        // RDKit❗✔️:   //     1 2
        // RDKit❗✔️:   //    /   \                                                   this space
        // RDKit❗✔️:   //   /     \                                       intentionally left blank
        // RDKit❗✔️:   //  a-------b
        // RDKit❗✔️:   //
        // RDKit❗✔️:   // If 1 and 2 are too close to each other we open up angle(1ab) if 1 is a
        // RDKit❗✔️:   // degree 1 node and
        // RDKit❗✔️:   // angle(2ba) if 2 is a degree 1 node. Say 1 is a degree 1 node but 2 is not.
        // RDKit❗✔️:   // Then from the neighbors of 2 we need to choose which one should be b. Also
        // RDKit❗✔️:   // keep in mind
        // RDKit❗✔️:   // that a need not be a neighbor of b. In this case we will pick b to be the
        // RDKit❗✔️:   // closest neighbor of a
        // RDKit❗✔️:   PRECONDITION(dp_mol, "");
        // RDKit❗✔️:   PRECONDITION(dmat, "");
        // RDKit❗✔️:   auto deg1 = getDepictDegree(dp_mol->getAtomWithIdx(aid1));
        // RDKit❗✔️:   auto deg2 = getDepictDegree(dp_mol->getAtomWithIdx(aid2));
        // RDKit❗✔️:   auto fixed1 = d_eatoms.at(aid1).df_fixed;
        // RDKit❗✔️:   auto fixed2 = d_eatoms.at(aid2).df_fixed;
        // RDKit❗✔️:   if ((deg1 > 1 || fixed1) && (deg2 > 1 || fixed2)) {
        // RDKit❗✔️:     return;
        // RDKit❗✔️:   }
        // RDKit❗✔️:   unsigned int aidA;
        // RDKit❗✔️:   unsigned int aidB;
        // RDKit❗✔️:   int type = 0;
        // RDKit❗✔️:   if ((deg1 == 1 && !fixed1) && (deg2 == 1 && !fixed2)) {
        // RDKit❗✔️:     aidA = _findDeg1Neighbor(dp_mol, aid1);
        // RDKit❗✔️:     aidB = _findDeg1Neighbor(dp_mol, aid2);
        // RDKit❗✔️:     type = 1;
        // RDKit❗✔️:   } else if ((deg1 == 1 && !fixed1) && (deg2 > 1 || fixed2)) {
        // RDKit❗✔️:     aidA = _findDeg1Neighbor(dp_mol, aid1);
        // RDKit❗✔️:     aidB = _findClosestNeighbor(dp_mol, dmat, aidA, aid2);
        // RDKit❗✔️:     type = 2;
        // RDKit❗✔️:   } else {
        // RDKit❗✔️:     aidB = _findDeg1Neighbor(dp_mol, aid2);
        // RDKit❗✔️:     aidA = _findClosestNeighbor(dp_mol, dmat, aidB, aid1);
        // RDKit❗✔️:     type = 3;
        // RDKit❗✔️:   }
        // RDKit❗✔️:   auto v2 = d_eatoms.at(aid1).loc - d_eatoms.at(aidA).loc;
        // RDKit❗✔️:   auto v1 = d_eatoms.at(aidB).loc - d_eatoms.at(aidA).loc;
        // RDKit❗✔️:   auto cross = (v1.x) * (v2.y) - (v1.y) * (v2.x);
        // RDKit❗✔️:   double angle;
        // RDKit❗✔️:   RDGeom::Transform2D trans1, trans2;
        // RDKit❗✔️:   switch (type) {
        // RDKit❗✔️:     case 1:
        // RDKit❗✔️:       angle = ANGLE_OPEN;
        // RDKit❗✔️:       if (cross < 0) {
        // RDKit❗✔️:         angle *= -1.0;
        // RDKit❗✔️:       }
        // RDKit❗✔️:       trans1.SetTransform(d_eatoms[aidA].loc, angle);
        // RDKit❗✔️:       trans2.SetTransform(d_eatoms[aidB].loc, -1.0 * angle);
        // RDKit❗✔️:       trans1.TransformPoint(d_eatoms[aid1].loc);
        // RDKit❗✔️:       trans2.TransformPoint(d_eatoms[aid2].loc);
        // RDKit❗✔️:       break;
        // RDKit❗✔️:     case 2:
        // RDKit❗✔️:       angle = 2.0 * ANGLE_OPEN;
        // RDKit❗✔️:       if (cross < 0) {
        // RDKit❗✔️:         angle *= -1.0;
        // RDKit❗✔️:       }
        // RDKit❗✔️:       trans1.SetTransform(d_eatoms[aidA].loc, angle);
        // RDKit❗✔️:       trans1.TransformPoint(d_eatoms[aid1].loc);
        // RDKit❗✔️:       break;
        // RDKit❗✔️:     case 3:
        // RDKit❗✔️:       angle = -2.0 * ANGLE_OPEN;
        // RDKit❗✔️:       if (cross < 0) {
        // RDKit❗✔️:         angle *= -1.0;
        // RDKit❗✔️:       }
        // RDKit❗✔️:       trans2.SetTransform(d_eatoms[aidB].loc, angle);
        // RDKit❗✔️:       trans2.TransformPoint(d_eatoms[aid2].loc);
        // RDKit❗✔️:       break;
        // RDKit❗✔️:     default:
        // RDKit❗✔️:       break;
        // RDKit❗✔️:   }
        // RDKit❗✔️: }
        // Behavior: preserve source three-way branch and signed 0.1222 radians.
        // Complexity: O(degree) neighbor lookup, two O(1) transforms at most.
        let degree1 = self.topology.adjacency.neighbors_of(aid1).len();
        let degree2 = self.topology.adjacency.neighbors_of(aid2).len();
        let atom1 = self
            .atoms
            .get(&aid1)
            .ok_or(FragmentError::AtomNotEmbedded { atom: aid1 })?;
        let atom2 = self
            .atoms
            .get(&aid2)
            .ok_or(FragmentError::AtomNotEmbedded { atom: aid2 })?;
        let (fixed1, fixed2) = (atom1.fixed, atom2.fixed);
        if (degree1 > 1 || fixed1) && (degree2 > 1 || fixed2) {
            return Ok(());
        }
        let (a, b, case) = if degree1 == 1 && !fixed1 && degree2 == 1 && !fixed2 {
            (
                self.degree_one_neighbor(aid1)?,
                self.degree_one_neighbor(aid2)?,
                1,
            )
        } else if degree1 == 1 && !fixed1 && (degree2 > 1 || fixed2) {
            let a = self.degree_one_neighbor(aid1)?;
            (a, self.closest_neighbor(distance, a, aid2), 2)
        } else {
            let b = self.degree_one_neighbor(aid2)?;
            (self.closest_neighbor(distance, b, aid1), b, 3)
        };
        let a_loc = self
            .atoms
            .get(&a)
            .ok_or(FragmentError::AtomNotEmbedded { atom: a })?
            .loc;
        let b_loc = self
            .atoms
            .get(&b)
            .ok_or(FragmentError::AtomNotEmbedded { atom: b })?
            .loc;
        let atom1_loc = self.atoms[&aid1].loc;
        let v2 = [atom1_loc[0] - a_loc[0], atom1_loc[1] - a_loc[1]];
        let v1 = [b_loc[0] - a_loc[0], b_loc[1] - a_loc[1]];
        let cross = v1[0] * v2[1] - v1[1] * v2[0];
        let sign = if cross < 0.0 { -1.0 } else { 1.0 };
        match case {
            1 => {
                let angle = 0.1222 * sign;
                self.atoms.get_mut(&aid1).expect("embedded").loc =
                    Transform2D::around(a_loc, angle).transform_point(atom1_loc);
                let atom2_loc = self.atoms[&aid2].loc;
                self.atoms.get_mut(&aid2).expect("embedded").loc =
                    Transform2D::around(b_loc, -angle).transform_point(atom2_loc);
            }
            2 => {
                self.atoms.get_mut(&aid1).expect("embedded").loc =
                    Transform2D::around(a_loc, 0.2444 * sign).transform_point(atom1_loc);
            }
            3 => {
                let atom2_loc = self.atoms[&aid2].loc;
                self.atoms.get_mut(&aid2).expect("embedded").loc =
                    Transform2D::around(b_loc, -0.2444 * sign).transform_point(atom2_loc);
            }
            _ => unreachable!(),
        }
        Ok(())
    }

    pub(crate) fn remove_collisions_open_angles(&mut self) -> Result<(), FragmentError> {
        // RDKit❗✔️: void EmbeddedFrag::removeCollisionsOpenAngles() {
        // RDKit❗✔️:   auto dmat = RDKit::MolOps::getDistanceMat(*dp_mol);
        // RDKit❗✔️:   // try opening up angles
        // RDKit❗✔️:   for (const auto &cpi : this->findCollisions(dmat, 0)) {
        // RDKit❗✔️:     // find out which of the two offending atoms we want to move
        // RDKit❗✔️:     // we will use the one with the smallest degree
        // RDKit❗✔️:     this->openAngles(dmat, cpi.first, cpi.second);
        // RDKit❗✔️:   }
        // RDKit❗✔️: }
        // Behavior: snapshot the atom-only collision list and process its order.
        // Complexity: one graph matrix and collision discovery; per pair O(degree).
        let distance = self.collision_distance_matrix()?;
        let collisions = self.find_collisions(&distance, false);
        for (first, second) in collisions {
            self.open_angles(&distance, first, second)?;
        }
        Ok(())
    }

    fn two_ring_neighbor_path(&self, aid: usize) -> Vec<(usize, [usize; 2])> {
        // RDKit❗✔️: void _recurseDegTwoRingAtoms(unsigned int aid, const RDKit::ROMol *mol,
        // RDKit❗✔️:                              RDKit::INT_VECT &rPath,
        // RDKit❗✔️:                              RDKit::INT_INT_VECT_MAP &nbrMap) {
        // RDKit❗✔️:   PRECONDITION(mol, "");
        // RDKit❗✔️:   // find all atoms along a path that have two ring atoms on them
        // RDKit❗✔️:   // aid is where will start looking and then we will recurse
        // RDKit❗✔️:   RDKit::INT_VECT nbrs;
        // RDKit❗✔️:   for (const auto bnd : mol->atomBonds(mol->getAtomWithIdx(aid))) {
        // RDKit❗✔️:     if (mol->getRingInfo()->numBondRings(bnd->getIdx())) {
        // RDKit❗✔️:       nbrs.push_back(bnd->getOtherAtomIdx(aid));
        // RDKit❗✔️:     }
        // RDKit❗✔️:   }
        // RDKit❗✔️:   if (nbrs.size() != 2) {
        // RDKit❗✔️:     return;
        // RDKit❗✔️:   } else {
        // RDKit❗✔️:     rPath.push_back(aid);
        // RDKit❗✔️:     nbrMap[aid] = nbrs;
        // RDKit❗✔️:     for (auto nbr : nbrs) {
        // RDKit❗✔️:       if (std::find(rPath.begin(), rPath.end(), nbr) == rPath.end()) {
        // RDKit❗✔️:         _recurseDegTwoRingAtoms(nbr, mol, rPath, nbrMap);
        // RDKit❗✔️:       }
        // RDKit❗✔️:     }
        // RDKit❗✔️:   }
        // RDKit❗✔️: }
        // Behavior: recurse in bond-table adjacency order only through atoms
        // having exactly two ring bonds; visit each accepted row once.
        // Complexity: source-shaped linear neighbor scan plus visited-vector
        // search per recursion, with no topology copy.
        fn visit(frag: &EmbeddedFrag<'_>, aid: usize, path: &mut Vec<(usize, [usize; 2])>) {
            let ring_neighbors = frag
                .topology
                .adjacency
                .neighbors_of(aid)
                .iter()
                .filter(|row| frag.rings.num_bond_rings(row.bond) != 0)
                .map(|row| row.atom_index)
                .collect::<Vec<_>>();
            if ring_neighbors.len() != 2 {
                return;
            }
            let neighbors = [ring_neighbors[0], ring_neighbors[1]];
            path.push((aid, neighbors));
            for neighbor in neighbors {
                if !path.iter().any(|(known, _)| *known == neighbor) {
                    visit(frag, neighbor, path);
                }
            }
        }
        let mut path = Vec::new();
        visit(self, aid, &mut path);
        path
    }

    fn count_non_ring_bonds(&self, aid: usize, path: &[AtomId]) -> usize {
        // RDKit❗✔️: unsigned int _anyNonRingBonds(unsigned int aid, RDKit::INT_LIST path,
        // RDKit❗✔️:                               const RDKit::ROMol *mol) {
        // RDKit❗✔️:   PRECONDITION(mol, "");
        // RDKit❗✔️:   // check if there are any non-ring bonds on the path starting at aid
        // RDKit❗✔️:   auto prev = aid;
        // RDKit❗✔️:   auto nOpen = 0u;
        // RDKit❗✔️:   for (auto pi : path) {
        // RDKit❗✔️:     const auto bond = mol->getBondBetweenAtoms(prev, pi);
        // RDKit❗✔️:     CHECK_INVARIANT(bond, "no bond found");
        // RDKit❗✔️:     if (!mol->getRingInfo()->numBondRings(bond->getIdx())) {
        // RDKit❗✔️:       ++nOpen;
        // RDKit❗✔️:     }
        // RDKit❗✔️:     prev = pi;
        // RDKit❗✔️:   }
        // RDKit❗✔️:   return nOpen;
        // RDKit❗✔️: }
        // Behavior: path excludes the starting row and counts all non-ring
        // edges, not only the first one.
        // Complexity: O(path length * adjacency degree), source equivalent.
        let mut previous = aid;
        let mut count = 0;
        for row in path {
            let next = row.index();
            let edge = self
                .topology
                .adjacency
                .neighbors_of(previous)
                .iter()
                .find(|edge| edge.atom_index == next)
                .expect("shortest path edge");
            if self.rings.num_bond_rings(edge.bond) == 0 {
                count += 1;
            }
            previous = next;
        }
        count
    }

    pub(crate) fn remove_collisions_shorten_bonds(&mut self) -> Result<(), FragmentError> {
        // RDKit❗✔️: void EmbeddedFrag::removeCollisionsShortenBonds() {
        // RDKit❗✔️:   auto dmat = RDKit::MolOps::getDistanceMat(*dp_mol);
        // RDKit❗✔️:   // if there are still some collision points left - flipping rotatable bonds
        // RDKit❗✔️:   // and opening angles is not doing it - we will try two last things
        // RDKit❗✔️:   //  - if all the bonds between the colliding atoms are rings bonds,
        // RDKit❗✔️:   //    we most likely have a collision within a bridged system (Issue 199).
        // RDKit❗✔️:   //    In this case we will try to find a path of colliding atoms (in one
        // RDKit❗✔️:   //    of the rings) and shorten all the bond in the path
        // RDKit❗✔️:   //  - on the other hand if we have non-ring bonds as well in the path
        // RDKit❗✔️:   //    between the colliding atoms we will simply shorten each one of
        // RDKit❗✔️:   //    them by a little bit.
        // RDKit❗✔️:   auto colls = this->findCollisions(dmat, 0);
        // RDKit❗✔️:   auto ncols = colls.size();
        // RDKit❗✔️:   auto iter = 0u;
        // RDKit❗✔️:   while (ncols && iter < MAX_COLL_ITERS) {
        // RDKit❗✔️:     const auto cAids = colls.front();
        // RDKit❗✔️:     // find out which of the two offending atoms we want to move
        // RDKit❗✔️:     // we will use the one with the smallest degree
        // RDKit❗✔️:     auto aid1 = cAids.first;
        // RDKit❗✔️:     auto aid2 = cAids.second;
        // RDKit❗✔️:     auto fixed1 = d_eatoms.at(aid1).df_fixed;
        // RDKit❗✔️:     auto fixed2 = d_eatoms.at(aid2).df_fixed;
        // RDKit❗✔️:     if (fixed1 && fixed2) {
        // RDKit❗✔️:       // both atoms are fixed, so there's nothing
        // RDKit❗✔️:       // we can do about this collision.
        // RDKit❗✔️:       colls.erase(colls.begin());
        // RDKit❗✔️:       ncols = colls.size();
        // RDKit❗✔️:       ++iter;
        // RDKit❗✔️:       continue;
        // RDKit❗✔️:     }
        // RDKit❗✔️:     auto deg1 = dp_mol->getAtomWithIdx(aid1)->getDegree();
        // RDKit❗✔️:     auto deg2 = dp_mol->getAtomWithIdx(aid2)->getDegree();
        // RDKit❗✔️:     if (fixed1 || (deg2 > deg1 && !fixed2)) {
        // RDKit❗✔️:       // reverse the order
        // RDKit❗✔️:       std::swap(deg1, deg2);
        // RDKit❗✔️:       std::swap(aid1, aid2);
        // RDKit❗✔️:       std::swap(fixed1, fixed2);
        // RDKit❗✔️:     }
        // RDKit❗✔️:     // now find the path between the two ends
        // RDKit❗✔️:     auto path = RDKit::MolOps::getShortestPath(*dp_mol, aid1, aid2);
        // RDKit❗✔️:     if (!path.size()) {
        // RDKit❗✔️:       // there's no path between the ends, so there's nothing
        // RDKit❗✔️:       // we can really do about this collision.
        // RDKit❗✔️:       colls.erase(colls.begin());
        // RDKit❗✔️:     } else {
        // RDKit❗✔️:       // aid1 is on the front of the path, pop it off:
        // RDKit❗✔️:       CHECK_INVARIANT(path.front() == aid1, "bad path head");
        // RDKit❗✔️:       path.pop_front();
        // RDKit❗✔️:       auto nOpen = _anyNonRingBonds(aid1, path, dp_mol);
        // RDKit❗✔️:       if (nOpen > 0) {
        // RDKit❗✔️:         if (deg1 == 1) {
        // RDKit❗✔️:           auto loc = d_eatoms.at(aid1).loc;
        // RDKit❗✔️:           auto aidA = _findDeg1Neighbor(dp_mol, aid1);
        // RDKit❗✔️:           loc -= d_eatoms[aidA].loc;
        // RDKit❗✔️:           loc *= .9;
        // RDKit❗✔️:           if (loc.length() > .75) {
        // RDKit❗✔️:             loc += d_eatoms[aidA].loc;
        // RDKit❗✔️:             d_eatoms[aid1].loc = loc;
        // RDKit❗✔️:           }
        // RDKit❗✔️:         }
        // RDKit❗✔️:         if (deg2 == 1 && !fixed2) {
        // RDKit❗✔️:           auto loc = d_eatoms.at(aid2).loc;
        // RDKit❗✔️:           auto aidA = _findDeg1Neighbor(dp_mol, aid2);
        // RDKit❗✔️:           loc -= d_eatoms[aidA].loc;
        // RDKit❗✔️:           loc *= .9;
        // RDKit❗✔️:           if (loc.length() > .75) {
        // RDKit❗✔️:             loc += d_eatoms[aidA].loc;
        // RDKit❗✔️:             d_eatoms[aid2].loc = loc;
        // RDKit❗✔️:           }
        // RDKit❗✔️:         }
        // RDKit❗✔️:       } else {
        // RDKit❗✔️:         // we probably have a bridged system
        // RDKit❗✔️:         // lets hope that aids has only two ring bond on it
        // RDKit❗✔️:         RDKit::INT_VECT rPath;
        // RDKit❗✔️:         RDKit::INT_INT_VECT_MAP nbrMap;
        // RDKit❗✔️:         _recurseDegTwoRingAtoms(aid1, dp_mol, rPath, nbrMap);
        // RDKit❗✔️:         if (rPath.size() == 0) {
        // RDKit❗✔️:           _recurseDegTwoRingAtoms(aid2, dp_mol, rPath, nbrMap);
        // RDKit❗✔️:         }
        // RDKit❗✔️:         // now we will take each of the atoms in rPath and
        // RDKit❗✔️:         // "move them in" a little bit this is what "move them
        // RDKit❗✔️:         //  in" means (what we need is hand drawn picture in the comments)
        // RDKit❗✔️:         // - let r1 and r2 be the ring neighbor of the current atom r0
        // RDKit❗✔️:         // - we will find the vector that bisects angle(r1, r0, r2)
        // RDKit❗✔️:         // - we will move r0 along this vector
        // RDKit❗✔️:         RDGeom::INT_POINT2D_MAP moveMap;
        // RDKit❗✔️:         for (auto rpi : rPath) {
        // RDKit❗✔️:           if (d_eatoms.at(rpi).df_fixed) {
        // RDKit❗✔️:             continue;
        // RDKit❗✔️:           }
        // RDKit❗✔️:           auto mv = d_eatoms[nbrMap[rpi][0]].loc;
        // RDKit❗✔️:           mv += d_eatoms[nbrMap.at(rpi)[1]].loc;
        // RDKit❗✔️:           mv *= 0.5;
        // RDKit❗✔️:           mv -= d_eatoms.at(rpi).loc;
        // RDKit❗✔️:           mv.normalize();
        // RDKit❗✔️:           mv *= COLLISION_THRES;
        // RDKit❗✔️:           moveMap[rpi] = mv;
        // RDKit❗✔️:         }
        // RDKit❗✔️:         for (auto rpi : rPath) {
        // RDKit❗✔️:           d_eatoms[rpi].loc += moveMap[rpi];
        // RDKit❗✔️:         }
        // RDKit❗✔️:       }
        // RDKit❗✔️:       colls = this->findCollisions(dmat, 0);
        // RDKit❗✔️:     }
        // RDKit❗✔️:     ncols = colls.size();
        // RDKit❗✔️:     ++iter;
        // RDKit❗✔️:   }
        // RDKit❗✔️: }
        // Behavior: source first-pair and 15-iteration cap, fixed/degree swap,
        // disconnected erase, terminal >0.75 contraction and staged ring moves.
        // Complexity: one core distance matrix and BFS per iteration; ring
        // path recursion and collision recomputation match source shape.
        let distance = self.collision_distance_matrix()?;
        let mut collisions = self.find_collisions(&distance, false);
        let mut iteration = 0;
        while !collisions.is_empty() && iteration < 15 {
            let (mut aid1, mut aid2) = collisions[0];
            let (mut fixed1, mut fixed2) = (self.atoms[&aid1].fixed, self.atoms[&aid2].fixed);
            if fixed1 && fixed2 {
                collisions.remove(0);
                iteration += 1;
                continue;
            }
            let (mut degree1, mut degree2) = (
                self.topology.adjacency.neighbors_of(aid1).len(),
                self.topology.adjacency.neighbors_of(aid2).len(),
            );
            if fixed1 || (degree2 > degree1 && !fixed2) {
                std::mem::swap(&mut degree1, &mut degree2);
                std::mem::swap(&mut aid1, &mut aid2);
                std::mem::swap(&mut fixed1, &mut fixed2);
            }
            let path = shortest_path(self.topology, AtomId::new(aid1), AtomId::new(aid2))?;
            if path.is_empty() {
                collisions.remove(0);
            } else {
                let non_ring = self.count_non_ring_bonds(aid1, &path[1..]);
                if non_ring > 0 {
                    if degree1 == 1 {
                        self.shorten_terminal_bond(aid1)?;
                    }
                    if degree2 == 1 && !fixed2 {
                        self.shorten_terminal_bond(aid2)?;
                    }
                } else {
                    let mut ring_path = self.two_ring_neighbor_path(aid1);
                    if ring_path.is_empty() {
                        ring_path = self.two_ring_neighbor_path(aid2);
                    }
                    let mut moves = BTreeMap::<usize, Point2>::new();
                    for &(aid, [n1, n2]) in &ring_path {
                        if self.atoms[&aid].fixed {
                            continue;
                        }
                        let first = self
                            .atoms
                            .get(&n1)
                            .ok_or(FragmentError::AtomNotEmbedded { atom: n1 })?
                            .loc;
                        let second = self
                            .atoms
                            .get(&n2)
                            .ok_or(FragmentError::AtomNotEmbedded { atom: n2 })?
                            .loc;
                        let center = self.atoms[&aid].loc;
                        let vector = [
                            (first[0] + second[0]) * 0.5 - center[0],
                            (first[1] + second[1]) * 0.5 - center[1],
                        ];
                        let length = vector[0].hypot(vector[1]);
                        // RDKit❗✔️:   void normalize() override {
                        // RDKit❗✔️:     double ln = this->length();
                        // RDKit❗✔️:     if (ln < zero_tolerance) {
                        // RDKit❗✔️:       throw std::runtime_error("Cannot normalize a zero length vector");
                        // RDKit❗✔️:     }
                        // RDKit❗✔️:     x /= ln;
                        // RDKit❗✔️:     y /= ln;
                        // RDKit❗✔️:   }
                        // Behavior: preserve its 1e-16 zero-tolerance error.
                        // Complexity: one norm and two divisions.
                        if length < 1.0e-16 {
                            return Err(FragmentError::CoincidentPoints);
                        }
                        moves.insert(aid, [vector[0] / length * 0.70, vector[1] / length * 0.70]);
                    }
                    for (aid, _) in ring_path {
                        let movement = moves.get(&aid).copied().unwrap_or([0.0, 0.0]);
                        let row = &mut self
                            .atoms
                            .get_mut(&aid)
                            .expect("ring-path atom embedded")
                            .loc;
                        row[0] += movement[0];
                        row[1] += movement[1];
                    }
                }
                collisions = self.find_collisions(&distance, false);
            }
            iteration += 1;
        }
        Ok(())
    }

    fn shorten_terminal_bond(&mut self, aid: usize) -> Result<(), FragmentError> {
        // RDKit❗✔️:           auto loc = d_eatoms.at(aid1).loc;
        // RDKit❗✔️:           auto aidA = _findDeg1Neighbor(dp_mol, aid1);
        // RDKit❗✔️:           loc -= d_eatoms[aidA].loc;
        // RDKit❗✔️:           loc *= .9;
        // RDKit❗✔️:           if (loc.length() > .75) {
        // RDKit❗✔️:             loc += d_eatoms[aidA].loc;
        // RDKit❗✔️:             d_eatoms[aid1].loc = loc;
        // RDKit❗✔️:           }
        // Behavior: strict >0.75 check is on the contracted displacement.
        // Complexity: one direct adjacency and coordinate lookup.
        let neighbor = self.degree_one_neighbor(aid)?;
        let parent = self
            .atoms
            .get(&neighbor)
            .ok_or(FragmentError::AtomNotEmbedded { atom: neighbor })?
            .loc;
        let current = self.atoms[&aid].loc;
        let vector = [
            (current[0] - parent[0]) * 0.9,
            (current[1] - parent[1]) * 0.9,
        ];
        if vector[0].hypot(vector[1]) > 0.75 {
            self.atoms.get_mut(&aid).expect("embedded").loc =
                [parent[0] + vector[0], parent[1] + vector[1]];
        }
        Ok(())
    }

    pub(crate) fn find_collisions(
        &mut self,
        distance: &DenseMatrix,
        include_bonds: bool,
    ) -> Vec<(usize, usize)> {
        // RDKit❗✔️: std::vector<PAIR_I_I> EmbeddedFrag::findCollisions(const double *dmat,
        // RDKit❗✔️:                                                    bool includeBonds) {
        // RDKit❗✔️:   // find a pair of atoms that are too close to each other
        // RDKit❗✔️:   std::vector<PAIR_I_I> res;
        // RDKit❗✔️:   for (auto &d_eatom : d_eatoms) {
        // RDKit❗✔️:     d_eatom.second.d_density = 0.0;
        // RDKit❗✔️:   }
        // RDKit❗✔️:   auto colThres2 = COLLISION_THRES * COLLISION_THRES;
        // RDKit❗✔️:   // if we a re dealing with non carbon atoms we will increase the collision
        // RDKit❗✔️:   // threshold. This is because only hetero atoms are typically drawn in a
        // RDKit❗✔️:   // depiction.
        // RDKit❗✔️:   double atomTypeFactor1, atomTypeFactor2;
        // RDKit❗✔️:   for (auto efi = d_eatoms.begin(); efi != d_eatoms.end(); ++efi) {
        // RDKit❗✔️:     atomTypeFactor1 = 1.0;
        // RDKit❗✔️:     if (dp_mol->getAtomWithIdx(efi->first)->getAtomicNum() != 6) {
        // RDKit❗✔️:       atomTypeFactor1 = HETEROATOM_COLL_SCALE;
        // RDKit❗✔️:     }
        // RDKit❗✔️:     for (auto efj = d_eatoms.begin(); efj != efi; ++efj) {
        // RDKit❗✔️:       atomTypeFactor2 = 1.0;
        // RDKit❗✔️:       if (dp_mol->getAtomWithIdx(efj->first)->getAtomicNum() != 6) {
        // RDKit❗✔️:         atomTypeFactor2 = HETEROATOM_COLL_SCALE;
        // RDKit❗✔️:       }
        // RDKit❗✔️:       auto ptj = efj->second.loc;
        // RDKit❗✔️:       ptj -= efi->second.loc;
        // RDKit❗✔️:       auto d2 = ptj.lengthSq();
        // RDKit❗✔️:       if (d2 > 1.0e-3) {
        // RDKit❗✔️:         efi->second.d_density += (1 / d2);
        // RDKit❗✔️:         efj->second.d_density += (1 / d2);
        // RDKit❗✔️:       } else {
        // RDKit❗✔️:         efi->second.d_density += 1000.0;
        // RDKit❗✔️:         efj->second.d_density += 1000.0;
        // RDKit❗✔️:       }
        // RDKit❗✔️:       d2 /= (atomTypeFactor1 * atomTypeFactor2);
        // RDKit❗✔️:       if (d2 < colThres2) {
        // RDKit❗✔️:         PAIR_I_I cAids(efi->first, efj->first);
        // RDKit❗✔️:         res.push_back(cAids);
        // RDKit❗✔️:       }
        // RDKit❗✔️:     }
        // RDKit❗✔️:   }
        // RDKit❗✔️:   if (includeBonds) {
        // RDKit❗✔️:     // now find bond collisions
        // RDKit❗✔️:     double BOND_THRES2 = BOND_THRES * BOND_THRES;
        // RDKit❗✔️:     for (const auto b1 : dp_mol->bonds()) {
        // RDKit❗✔️:       auto bid1 = b1->getIdx();
        // RDKit❗✔️:       auto beg1 = b1->getBeginAtomIdx();
        // RDKit❗✔️:       auto end1 = b1->getEndAtomIdx();
        // RDKit❗✔️:       if ((d_eatoms.find(beg1) != d_eatoms.end()) &&
        // RDKit❗✔️:           (d_eatoms.find(end1) != d_eatoms.end())) {
        // RDKit❗✔️:         auto v1 = d_eatoms[end1].loc - d_eatoms[beg1].loc;
        // RDKit❗✔️:         auto avg1 = d_eatoms[end1].loc + d_eatoms[beg1].loc;
        // RDKit❗✔️:         avg1 *= 0.5;
        // RDKit❗✔️:         for (const auto b2 : dp_mol->bonds()) {
        // RDKit❗✔️:           if (b2->getIdx() <= bid1) {
        // RDKit❗✔️:             continue;
        // RDKit❗✔️:           }
        //
        // RDKit❗✔️:           auto beg2 = b2->getBeginAtomIdx();
        // RDKit❗✔️:           auto end2 = b2->getEndAtomIdx();
        // RDKit❗✔️:           if ((d_eatoms.find(beg2) != d_eatoms.end()) &&
        // RDKit❗✔️:               (d_eatoms.find(end2) != d_eatoms.end())) {
        // RDKit❗✔️:             auto avg2 = d_eatoms[end2].loc + d_eatoms[beg2].loc;
        // RDKit❗✔️:             avg2 *= 0.5;
        // RDKit❗✔️:             avg2 -= avg1;
        // RDKit❗✔️:             if (avg2.lengthSq() < 0.5 && avg2.lengthSq() < BOND_THRES2) {
        // RDKit❗✔️:               auto v2 = d_eatoms[beg2].loc - d_eatoms[beg1].loc;
        // RDKit❗✔️:               auto v3 = d_eatoms[end2].loc - d_eatoms[beg1].loc;
        // RDKit❗✔️:               auto valProd = _crossVal(v1, v2) * _crossVal(v1, v3);
        // RDKit❗✔️:               if (valProd < -1e-6) {
        // RDKit❗✔️:                 // we have a collision, find the closest two atoms
        // RDKit❗✔️:                 auto cAids =
        // RDKit❗✔️:                     _findClosestPair(beg1, end1, beg2, end2, *dp_mol, dmat);
        // RDKit❗✔️:                 res.push_back(cAids);
        // RDKit❗✔️:               }
        // RDKit❗✔️:             }
        // RDKit❗✔️:           }
        // RDKit❗✔️:         }
        // RDKit❗✔️:       }
        // RDKit❗✔️:     }
        // RDKit❗✔️:   }
        // RDKit❗✔️:   return res;
        // RDKit❗✔️: }
        // RDKit❗✔️: PAIR_I_I _findClosestPair(unsigned int beg1, unsigned int end1,
        // RDKit❗✔️:                           unsigned int beg2, unsigned int end2,
        // RDKit❗✔️:                           const RDKit::ROMol &mol, const double *dmat) {
        // RDKit❗✔️:   auto na = mol.getNumAtoms();
        // RDKit❗✔️:   auto d1 = dmat[beg1 * na + beg2];
        // RDKit❗✔️:   auto d2 = dmat[beg1 * na + end2];
        // RDKit❗✔️:   auto d3 = dmat[end1 * na + beg2];
        // RDKit❗✔️:   auto d4 = dmat[end1 * na + end2];
        // RDKit❗✔️:   auto minPr =
        // RDKit❗✔️:       std::min(PAIR_D_I_I(d1, PAIR_I_I(beg1, beg2)),
        // RDKit❗✔️:                PAIR_D_I_I(d2, PAIR_I_I(beg1, end2)), _pairDIICompAscending);
        // RDKit❗✔️:   minPr = std::min(minPr, PAIR_D_I_I(d3, PAIR_I_I(end1, beg2)),
        // RDKit❗✔️:                    _pairDIICompAscending);
        // RDKit❗✔️:   minPr = std::min(minPr, PAIR_D_I_I(d4, PAIR_I_I(end1, end2)),
        // RDKit❗✔️:                    _pairDIICompAscending);
        // RDKit❗✔️:   return minPr.second;
        // RDKit❗✔️: }
        // RDKit❗✔️: int _pairDIICompAscending(const PAIR_D_I_I &arg1, const PAIR_D_I_I &arg2) {
        // RDKit❗✔️:   return (arg1.first < arg2.first);
        // RDKit❗✔️: }
        // RDKit❗✔️: double _crossVal(const RDGeom::Point2D &v1, const RDGeom::Point2D &v2) {
        // RDKit❗✔️:   return v1.x * v2.y - v2.x * v1.y;
        // RDKit❗✔️: }
        // Behavior: source BTreeMap key order, strict thresholds, four-pair
        // closest comparison and bond-table order are retained. Query carrier
        // elements are not replaced by inferred predicate atoms.
        // Complexity: O(embedded atoms^2 + topology bonds^2), as source;
        // no topology or coordinate-block clone.
        for atom in self.atoms.values_mut() {
            atom.density = 0.0;
        }
        let keys: Vec<_> = self.atoms.keys().copied().collect();
        let mut collisions = Vec::new();
        for (outer, &aid) in keys.iter().enumerate() {
            let factor1 = if self.topology.atoms[aid].atomic_number() != 6 {
                1.3
            } else {
                1.0
            };
            for &bid in &keys[..outer] {
                let factor2 = if self.topology.atoms[bid].atomic_number() != 6 {
                    1.3
                } else {
                    1.0
                };
                let a = self.atoms[&aid].loc;
                let b = self.atoms[&bid].loc;
                let d2 = (a[0] - b[0]).powi(2) + (a[1] - b[1]).powi(2);
                let increment = if d2 > 1.0e-3 { 1.0 / d2 } else { 1000.0 };
                self.atoms.get_mut(&aid).expect("embedded").density += increment;
                self.atoms.get_mut(&bid).expect("embedded").density += increment;
                if d2 / (factor1 * factor2) < 0.70 * 0.70 {
                    collisions.push((aid, bid));
                }
            }
        }
        if include_bonds {
            for (first_index, first) in self.topology.bonds.iter().enumerate() {
                let (a, b) = (first.begin().index(), first.end().index());
                let (Some(pa), Some(pb)) = (self.atoms.get(&a), self.atoms.get(&b)) else {
                    continue;
                };
                let v1 = [pb.loc[0] - pa.loc[0], pb.loc[1] - pa.loc[1]];
                let midpoint1 = [(pa.loc[0] + pb.loc[0]) * 0.5, (pa.loc[1] + pb.loc[1]) * 0.5];
                for second in &self.topology.bonds[first_index + 1..] {
                    let (c, d) = (second.begin().index(), second.end().index());
                    let (Some(pc), Some(pd)) = (self.atoms.get(&c), self.atoms.get(&d)) else {
                        continue;
                    };
                    let midpoint2 = [(pc.loc[0] + pd.loc[0]) * 0.5, (pc.loc[1] + pd.loc[1]) * 0.5];
                    let midpoint_distance2 = (midpoint2[0] - midpoint1[0]).powi(2)
                        + (midpoint2[1] - midpoint1[1]).powi(2);
                    if midpoint_distance2 < 0.5 && midpoint_distance2 < 0.50 * 0.50 {
                        let v2 = [pc.loc[0] - pa.loc[0], pc.loc[1] - pa.loc[1]];
                        let v3 = [pd.loc[0] - pa.loc[0], pd.loc[1] - pa.loc[1]];
                        let cross2 = v1[0] * v2[1] - v2[0] * v1[1];
                        let cross3 = v1[0] * v3[1] - v3[0] * v1[1];
                        if cross2 * cross3 < -1.0e-6 {
                            let candidates = [(a, c), (a, d), (b, c), (b, d)];
                            let mut chosen = candidates[0];
                            let mut best = distance.get(chosen.0, chosen.1).expect("full matrix");
                            for candidate in &candidates[1..] {
                                let value =
                                    distance.get(candidate.0, candidate.1).expect("full matrix");
                                if value < best {
                                    chosen = *candidate;
                                    best = value;
                                }
                            }
                            collisions.push(chosen);
                        }
                    }
                }
            }
        }
        collisions
    }

    pub(crate) fn total_density(&self) -> f64 {
        // RDKit❗✔️: double EmbeddedFrag::totalDensity() {
        // RDKit❗✔️:   return std::accumulate(
        // RDKit❗✔️:       d_eatoms.begin(), d_eatoms.end(), 0.0,
        // RDKit❗✔️:       [](double accum, auto &dea) { return dea.second.d_density + accum; });
        // RDKit❗✔️: }
        // Behavior: source key-order left fold preserves floating summation.
        // Complexity: one scan, no allocation.
        self.atoms
            .values()
            .fold(0.0, |sum, atom| atom.density + sum)
    }

    pub(crate) fn rotatable_bonds_on_shortest_path(
        &self,
        aid1: usize,
        aid2: usize,
    ) -> Result<Vec<usize>, FragmentError> {
        // RDKit❗✔️: RDKit::INT_VECT getRotatableBonds(const RDKit::ROMol &mol, unsigned int aid1,
        // RDKit❗✔️:                                   unsigned int aid2) {
        // RDKit❗✔️:   PRECONDITION(aid1 < mol.getNumAtoms(), "");
        // RDKit❗✔️:   PRECONDITION(aid2 < mol.getNumAtoms(), "");
        // RDKit❗✔️:   RDKit::INT_LIST path = RDKit::MolOps::getShortestPath(mol, aid1, aid2);
        // RDKit❗✔️:   RDKit::INT_VECT res;
        // RDKit❗✔️:   if (path.size() >= 4) {
        // RDKit❗✔️:     // remove the first atom (aid1) and last atom (aid2)
        // RDKit❗✔️:     CHECK_INVARIANT(static_cast<unsigned int>(path.front()) == aid1,
        // RDKit❗✔️:                     "bad first element");
        // RDKit❗✔️:     path.pop_front();
        // RDKit❗✔️:     CHECK_INVARIANT(static_cast<unsigned int>(path.back()) == aid2,
        // RDKit❗✔️:                     "bad last element");
        // RDKit❗✔️:     path.pop_back();
        // RDKit❗✔️:     auto pid = path.front();
        // RDKit❗✔️:     for (auto aid : path) {
        // RDKit❗✔️:       if (aid == pid) {
        // RDKit❗✔️:         continue;
        // RDKit❗✔️:       }
        // RDKit❗✔️:       const RDKit::Bond *bond = mol.getBondBetweenAtoms(pid, aid);
        // RDKit❗✔️:       int bid = bond->getIdx();
        // RDKit❗✔️:       if ((bond->getStereo() <= RDKit::Bond::STEREOANY) &&
        // RDKit❗✔️:           (!(mol.getRingInfo()->numBondRings(bid)))) {
        // RDKit❗✔️:         res.push_back(bid);
        // RDKit❗✔️:       }
        // RDKit❗✔️:       pid = aid;
        // RDKit❗✔️:     }
        // RDKit❗✔️:   }
        // RDKit❗✔️:   return res;
        // RDKit❗✔️: }
        // Behavior: canonical core shortest path supplies source BFS ties;
        // only its interior edges are candidates, with no bond-order filter.
        // Complexity: core BFS plus O(path length * adjacency degree).
        let path = shortest_path(self.topology, AtomId::new(aid1), AtomId::new(aid2))?;
        let mut result = Vec::new();
        if path.len() >= 4 {
            for pair in path[1..path.len() - 1].windows(2) {
                let first = pair[0].index();
                let next = pair[1].index();
                let neighbor = self
                    .topology
                    .adjacency
                    .neighbors_of(first)
                    .iter()
                    .find(|row| row.atom_index == next)
                    .expect("shortest path edge");
                let bond = &self.topology.bonds[neighbor.bond.index()];
                if bond.stereo().rdkit_code() <= BondStereo::Any.rdkit_code()
                    && self.rings.num_bond_rings(bond.id()) == 0
                {
                    result.push(bond.id().index());
                }
            }
        }
        Ok(result)
    }
}

pub(crate) fn embed_cis_trans_systems<'a>(
    topology: &'a TopologyBlock,
    rings: &'a RingInfo,
) -> Result<Vec<EmbeddedFrag<'a>>, FragmentError> {
    // RDKit❗✔️: void embedCisTransSystems(const RDKit::ROMol &mol,
    // RDKit❗✔️:                           std::list<EmbeddedFrag> &efrags) {
    // RDKit❗✔️:   for (auto bond : mol.bonds()) {
    // RDKit❗✔️:     // check if this bond is in a cis/trans double bond
    // RDKit❗✔️:     // and it is not a ring bond
    // RDKit❗✔️:     if ((bond->getBondType() == RDKit::Bond::DOUBLE)  // this is a double bond
    // RDKit❗✔️:         && (bond->getStereo() >
    // RDKit❗✔️:             RDKit::Bond::STEREOANY)  // and has stereo chemistry specified
    // RDKit❗✔️:         && (!bond->getOwningMol().getRingInfo()->numBondRings(
    // RDKit❗✔️:                bond->getIdx()))) {  // not in a ring
    // RDKit❗✔️:       if (bond->getStereoAtoms().size() != 2) {
    // RDKit❗✔️:         BOOST_LOG(rdWarningLog)
    // RDKit❗✔️:             << "WARNING: bond found with stereo spec but no stereo atoms"
    // RDKit❗✔️:             << std::endl;
    // RDKit❗✔️:         continue;
    // RDKit❗✔️:       }
    // RDKit❗✔️:       EmbeddedFrag efrag(bond);
    // RDKit❗✔️:       efrag.setupNewNeighs();
    // RDKit❗✔️:       efrags.push_back(efrag);
    // RDKit❗✔️:     }
    // RDKit❗✔️:   }
    // RDKit❗✔️: }
    // Behavior: source bond-table order and missing-stereo-row skip.
    // Complexity: one bond scan, each selected seed updates its two neighbors.
    let mut fragments = Vec::new();
    for bond in &topology.bonds {
        if bond.order() != BondOrder::Double
            || bond.stereo().rdkit_code() <= BondStereo::Any.rdkit_code()
            || rings.num_bond_rings(bond.id()) != 0
            || bond.stereo_atoms().is_none()
        {
            continue;
        }
        let mut fragment = EmbeddedFrag::from_cis_trans_bond(bond.id(), topology, rings)?;
        fragment.setup_new_neighbors()?;
        fragments.push(fragment);
    }
    Ok(fragments)
}

pub(crate) fn embed_fused_systems<'a>(
    topology: &'a TopologyBlock,
    rings: &'a RingInfo,
    coordinate_map: Option<&PointMap>,
    use_templates: bool,
    templates: &mut CoordinateTemplates,
) -> Result<Vec<EmbeddedFrag<'a>>, FragmentError> {
    // RDKit❗❗: void embedFusedSystems(const RDKit::ROMol &mol,
    // RDKit❗❗:                        const RDKit::VECT_INT_VECT &arings,
    // RDKit❗❗:                        std::list<EmbeddedFrag> &efrags,
    // RDKit❗❗:                        const RDGeom::INT_POINT2D_MAP *coordMap,
    // RDKit❗❗:                        bool useRingTemplates) {
    // RDKit❗❗:   RDKit::INT_INT_VECT_MAP neighMap;
    // RDKit❗❗:   RingUtils::makeRingNeighborMap(arings, neighMap);
    // RDKit❗❗:   auto cnrs = arings.size();
    // RDKit❗❗:   boost::dynamic_bitset<> fusDone(cnrs);
    // RDKit❗❗:   auto curr = 0u;
    // RDKit❗❗:   while (curr < cnrs) {
    // RDKit❗❗:     RDKit::INT_VECT fused;
    // RDKit❗❗:     RingUtils::pickFusedRings(curr, neighMap, fused, fusDone);
    // RDKit❗❗:     RDKit::VECT_INT_VECT frings;
    // RDKit❗❗:     frings.reserve(fused.size());
    // RDKit❗❗:     for (auto rid : fused) {
    // RDKit❗❗:       frings.push_back(arings.at(rid));
    // RDKit❗❗:     }
    // RDKit❗❗:     bool allowRingTemplates = useRingTemplates;
    // RDKit❗❗:     if (useRingTemplates && coordMap) {
    // RDKit❗❗:       boost::dynamic_bitset<> coordMapAtoms(mol.getNumAtoms());
    // RDKit❗❗:       for (const auto &ring : frings) {
    // RDKit❗❗:         for (const auto &aid : ring) {
    // RDKit❗❗:           if (coordMap->find(aid) != coordMap->end()) {
    // RDKit❗❗:             coordMapAtoms.set(aid);
    // RDKit❗❗:           }
    // RDKit❗❗:         }
    // RDKit❗❗:       }
    // RDKit❗❗:       allowRingTemplates = (coordMapAtoms.count() < 2);
    // RDKit❗❗:     }
    // RDKit❗❗:     EmbeddedFrag efrag(&mol, frings, allowRingTemplates);
    // RDKit❗❗:     efrag.setupNewNeighs();
    // RDKit❗❗:     efrags.push_back(efrag);
    // RDKit❗❗:     size_t rix;
    // RDKit❗❗:     for (rix = 0; rix < cnrs; ++rix) {
    // RDKit❗❗:       if (!fusDone[rix]) {
    // RDKit❗❗:         curr = rix;
    // RDKit❗❗:         break;
    // RDKit❗❗:       }
    // RDKit❗❗:     }
    // RDKit❗❗:     if (rix == cnrs) {
    // RDKit❗❗:       break;
    // RDKit❗❗:     }
    // RDKit❗❗:   }
    // RDKit❗❗: }
    // Behavior: DFS source ring order and distinct constrained atoms; this
    // remains a private detached fragment producer, not final layout.
    // Complexity: one ring-neighbor graph and one fragment per system.
    let atom_rings: Vec<Vec<usize>> = rings
        .atom_rings()
        .iter()
        .map(|ring| ring.iter().map(|atom| atom.index()).collect())
        .collect();
    let mut fragments = Vec::new();
    for system in ring_systems(&atom_rings) {
        let fused_rings: Vec<_> = system
            .iter()
            .map(|&index| atom_rings[index].clone())
            .collect();
        let mut allow_templates = use_templates;
        if use_templates {
            if let Some(map) = coordinate_map {
                let constrained = ring_union(&fused_rings)
                    .iter()
                    .filter(|atom| map.contains_key(atom))
                    .count();
                allow_templates = constrained < 2;
            }
        }
        let mut fragment = EmbeddedFrag::from_fused_rings(
            topology,
            rings,
            &fused_rings,
            allow_templates,
            templates,
        )?;
        fragment.setup_new_neighbors()?;
        fragments.push(fragment);
    }
    Ok(fragments)
}

impl<'a> EmbeddedFrag<'a> {
    pub(crate) fn merge_frags_with_common(
        &mut self,
        fragments: &mut Vec<Self>,
    ) -> Result<(), FragmentError> {
        // RDKit❗✔️: void EmbeddedFrag::mergeFragsWithComm(std::list<EmbeddedFrag> &efrags) {
        // RDKit❗✔️:   PRECONDITION(dp_mol, "");
        // RDKit❗✔️:   auto nfri = efrags.end();
        // RDKit❗✔️:   while (1) {
        // RDKit❗✔️:     RDKit::INT_VECT commAtms;
        // RDKit❗✔️:     for (auto efri = efrags.begin(); efri != efrags.end(); ++efri) {
        // RDKit❗✔️:       if (!efri->isDone()) {
        // RDKit❗✔️:         commAtms = this->findCommonAtoms(*efri);
        // RDKit❗✔️:         if (commAtms.size() > 0) {
        // RDKit❗✔️:           nfri = efri;
        // RDKit❗✔️:           break;
        // RDKit❗✔️:         }
        // RDKit❗✔️:       }
        // RDKit❗✔️:     }
        // RDKit❗✔️:     if (commAtms.empty()) {
        // RDKit❗✔️:       break;
        // RDKit❗✔️:     }
        // RDKit❗✔️:     CHECK_INVARIANT(nfri != efrags.end(), "iterator not initialized");
        // RDKit❗✔️:     this->mergeWithCommon((*nfri), commAtms);
        // RDKit❗✔️:     for (auto cai : commAtms) {
        // RDKit❗✔️:       if (d_eatoms.at(cai).neighs.empty() &&
        // RDKit❗✔️:           (std::find(d_attachPts.begin(), d_attachPts.end(), cai) !=
        // RDKit❗✔️:            d_attachPts.end())) {
        // RDKit❗✔️:         d_attachPts.erase(
        // RDKit❗✔️:             std::remove(d_attachPts.begin(), d_attachPts.end(), cai));
        // RDKit❗✔️:       }
        // RDKit❗✔️:     }
        // RDKit❗✔️:     efrags.erase(nfri);
        // RDKit❗✔️:   }
        // RDKit❗✔️: }
        loop {
            let next = fragments
                .iter()
                .enumerate()
                .filter(|(_, fragment)| !fragment.done)
                .find_map(|(index, fragment)| {
                    let common = self.find_common_atoms(fragment);
                    (!common.is_empty()).then_some((index, common))
                });
            let Some((index, mut common)) = next else {
                break;
            };
            self.merge_with_common(&mut fragments[index], &mut common)?;
            for aid in common {
                if self
                    .atoms
                    .get(&aid)
                    .ok_or(FragmentError::AtomNotEmbedded { atom: aid })?
                    .neighs
                    .is_empty()
                {
                    self.attachment_points.retain(|&entry| entry != aid);
                }
            }
            fragments.remove(index);
        }
        Ok(())
    }

    pub(crate) fn expand_fragment(
        &mut self,
        nonring_atoms: &mut Vec<usize>,
        fragments: &mut Vec<Self>,
    ) -> Result<(), FragmentError> {
        // RDKit❗✔️: void EmbeddedFrag::expandEfrag(RDKit::INT_LIST &nratms,
        // RDKit❗✔️:                                std::list<EmbeddedFrag> &efrags) {
        // RDKit❗✔️:   PRECONDITION(dp_mol, "");
        // RDKit❗✔️:   this->mergeFragsWithComm(efrags);
        // RDKit❗✔️:   while (d_attachPts.size() > 0) {
        // RDKit❗✔️:     auto aid = d_attachPts.front();
        // RDKit❗✔️:     auto nbrs = d_eatoms[aid].neighs;
        // RDKit❗✔️:     CHECK_INVARIANT(!nbrs.empty(), "");
        // RDKit❗✔️:     for (auto nbri : nbrs) {
        // RDKit❗✔️:       auto nratmi = std::find(nratms.begin(), nratms.end(), nbri);
        // RDKit❗✔️:       if (nratmi != nratms.end()) {
        // RDKit❗✔️:         this->addNonRingAtom(nbri, aid);
        // RDKit❗✔️:         nratms.erase(nratmi);
        // RDKit❗✔️:       } else {
        // RDKit❗✔️:         auto nfri = efrags.end();
        // RDKit❗✔️:         for (auto efri = efrags.begin(); efri != efrags.end(); ++efri) {
        // RDKit❗✔️:           if (!efri->isDone()) {
        // RDKit❗✔️:             const auto &eatoms = efri->GetEmbeddedAtoms();
        // RDKit❗✔️:             if (eatoms.find(nbri) != eatoms.end()) {
        // RDKit❗✔️:               nfri = efri;
        // RDKit❗✔️:               break;
        // RDKit❗✔️:             }
        // RDKit❗✔️:           }
        // RDKit❗✔️:         }
        // RDKit❗✔️:         if (nfri != efrags.end()) {
        // RDKit❗✔️:           this->mergeNoCommon((*nfri), aid, nbri);
        // RDKit❗✔️:           if (d_eatoms.at(nbri).neighs.empty() &&
        // RDKit❗✔️:               (std::find(d_attachPts.begin(), d_attachPts.end(), nbri) !=
        // RDKit❗✔️:                d_attachPts.end())) {
        // RDKit❗✔️:             d_attachPts.erase(
        // RDKit❗✔️:                 std::remove(d_attachPts.begin(), d_attachPts.end(), nbri));
        // RDKit❗✔️:           }
        // RDKit❗✔️:           efrags.erase(nfri);
        // RDKit❗✔️:         }
        // RDKit❗✔️:       }
        // RDKit❗✔️:     }
        // RDKit❗✔️:     d_attachPts.pop_front();
        // RDKit❗✔️:     d_eatoms[aid].neighs.clear();
        // RDKit❗✔️:     this->mergeFragsWithComm(efrags);
        // RDKit❗✔️:   }
        // RDKit❗✔️: }
        self.merge_frags_with_common(fragments)?;
        while let Some(&aid) = self.attachment_points.first() {
            let neighbors = self
                .atoms
                .get(&aid)
                .ok_or(FragmentError::AtomNotEmbedded { atom: aid })?
                .neighs
                .clone();
            if neighbors.is_empty() {
                return Err(FragmentError::EmptyAttachment { atom: aid });
            }
            for neighbor in neighbors {
                if let Some(index) = nonring_atoms.iter().position(|&row| row == neighbor) {
                    self.add_non_ring_atom(neighbor, aid)?;
                    nonring_atoms.remove(index);
                } else if let Some(index) = fragments
                    .iter()
                    .position(|fragment| !fragment.done && fragment.atoms.contains_key(&neighbor))
                {
                    self.merge_no_common(&mut fragments[index], aid, neighbor)?;
                    if self
                        .atoms
                        .get(&neighbor)
                        .ok_or(FragmentError::AtomNotEmbedded { atom: neighbor })?
                        .neighs
                        .is_empty()
                    {
                        self.attachment_points.retain(|&entry| entry != neighbor);
                    }
                    fragments.remove(index);
                }
            }
            self.attachment_points.remove(0);
            self.atoms
                .get_mut(&aid)
                .ok_or(FragmentError::AtomNotEmbedded { atom: aid })?
                .neighs
                .clear();
            self.merge_frags_with_common(fragments)?;
        }
        Ok(())
    }

    pub(crate) fn find_common_atoms(&self, other: &Self) -> Vec<usize> {
        // RDKit❗✔️: RDKit::INT_VECT EmbeddedFrag::findCommonAtoms(const EmbeddedFrag &efrag2) {
        // RDKit❗✔️:   RDKit::INT_VECT res;
        // RDKit❗✔️:   for (auto eri1 : this->GetEmbeddedAtoms()) {
        // RDKit❗✔️:     for (auto eri2 : efrag2.GetEmbeddedAtoms()) {
        // RDKit❗✔️:       if (eri1.first == eri2.first) {
        // RDKit❗✔️:         res.push_back(eri1.first);
        // RDKit❗✔️:       }
        // RDKit❗✔️:     }
        // RDKit❗✔️:   }
        // RDKit❗✔️:   return res;
        // RDKit❗✔️: }
        let mut common = Vec::new();
        for aid in self.atoms.keys().copied() {
            for other_aid in other.atoms.keys().copied() {
                if aid == other_aid {
                    common.push(aid);
                }
            }
        }
        common
    }

    pub(crate) fn merge_no_common(
        &mut self,
        other: &mut Self,
        to_aid: usize,
        nbr_aid: usize,
    ) -> Result<(), FragmentError> {
        // RDKit❗✔️: void EmbeddedFrag::mergeNoCommon(EmbeddedFrag &embObj, unsigned int toAid,
        // RDKit❗✔️:                                  unsigned int nbrAid) {
        // RDKit❗✔️:   PRECONDITION(dp_mol, "");
        // RDKit❗✔️:   PRECONDITION(dp_mol == embObj.getMol(), "Molecule mismatch");
        // RDKit❗✔️:   RDKit::INT_VECT commAtms;
        // RDKit❗✔️:   this->addNonRingAtom(nbrAid, toAid);
        // RDKit❗✔️:   embObj.addNonRingAtom(toAid, nbrAid);
        // RDKit❗✔️:   commAtms.push_back(toAid);
        // RDKit❗✔️:   commAtms.push_back(nbrAid);
        // RDKit❗✔️:   this->mergeWithCommon(embObj, commAtms);
        // RDKit❗✔️: }
        self.check_same_topology(other)?;
        self.add_non_ring_atom(nbr_aid, to_aid)?;
        other.add_non_ring_atom(to_aid, nbr_aid)?;
        self.merge_with_common(other, &mut vec![to_aid, nbr_aid])
    }

    pub(crate) fn merge_with_common(
        &mut self,
        other: &mut Self,
        common: &mut Vec<usize>,
    ) -> Result<(), FragmentError> {
        // RDKit❗✔️: void EmbeddedFrag::mergeWithCommon(EmbeddedFrag &embObj,
        // RDKit❗✔️:                                    RDKit::INT_VECT &commAtms) {
        // RDKit❗✔️:   PRECONDITION(dp_mol, "");
        // RDKit❗✔️:   PRECONDITION(dp_mol == embObj.getMol(), "Molecule mismatch");
        // RDKit❗✔️:   PRECONDITION(commAtms.size() >= 1, "");
        // RDKit❗✔️:   unsigned int ctCase =
        // RDKit❗✔️:       0;
        // RDKit❗✔️:   if (commAtms.size() == 1) {
        // RDKit❗✔️:     auto commAid = commAtms.front();
        // RDKit❗✔️:     int otherAtom = -1;
        // RDKit❗✔️:     if (d_eatoms[commAid].CisTransNbr >= 0) {
        // RDKit❗✔️:       ctCase = 2;
        // RDKit❗✔️:       otherAtom =
        // RDKit❗✔️:           d_eatoms[commAid].nbr1;
        // RDKit❗✔️:       embObj.addNonRingAtom(otherAtom, commAid);
        // RDKit❗✔️:     } else if (embObj.d_eatoms[commAid].CisTransNbr >= 0) {
        // RDKit❗✔️:       ctCase = 1;
        // RDKit❗✔️:       otherAtom = embObj.d_eatoms[commAid].nbr1;
        // RDKit❗✔️:       this->addNonRingAtom(otherAtom, commAid);
        // RDKit❗✔️:     } else {
        // RDKit❗✔️:       otherAtom = d_eatoms[commAid].nbr1;
        // RDKit❗✔️:       if (otherAtom >= 0) {
        // RDKit❗✔️:         embObj.addNonRingAtom(otherAtom, commAid);
        // RDKit❗✔️:       }
        // RDKit❗✔️:     }
        // RDKit❗✔️:     if (otherAtom >= 0) {
        // RDKit❗✔️:       commAtms.push_back(otherAtom);
        // RDKit❗✔️:     }
        // RDKit❗✔️:   }
        // RDKit❗✔️:   RDGeom::Transform2D rtrans;
        // RDKit❗✔️:   if (commAtms.size() == 1) {
        // RDKit❗✔️:     rtrans.assign(this->computeOneAtomTrans(commAtms.front(), embObj));
        // RDKit❗✔️:   } else {
        // RDKit❗✔️:     auto cid1 = commAtms[0];
        // RDKit❗✔️:     auto cid2 = commAtms[1];
        // RDKit❗✔️:     const auto &ref1 = d_eatoms.at(cid1).loc;
        // RDKit❗✔️:     const auto &ref2 = d_eatoms.at(cid2).loc;
        // RDKit❗✔️:     const auto &oth1 = embObj.GetEmbeddedAtom(cid1).loc;
        // RDKit❗✔️:     const auto &oth2 = embObj.GetEmbeddedAtom(cid2).loc;
        // RDKit❗✔️:     rtrans.SetTransform(ref1, ref2, oth1, oth2);
        // RDKit❗✔️:   }
        // RDKit❗✔️:   embObj.Transform(rtrans);
        // RDKit❗✔️:   if (commAtms.size() >= 2) {
        // RDKit❗✔️:     if (ctCase > 0) {
        // RDKit❗✔️:       reflectIfNecessaryCisTrans(embObj, ctCase, commAtms[0], commAtms[1]);
        // RDKit❗✔️:     } else if (commAtms.size() == 2) {
        // RDKit❗✔️:       reflectIfNecessaryDensity(embObj, commAtms[0], commAtms[1]);
        // RDKit❗✔️:     } else {
        // RDKit❗✔️:       reflectIfNecessaryThirdPt(embObj, commAtms[0], commAtms[1], commAtms[2]);
        // RDKit❗✔️:     }
        // RDKit❗✔️:   }
        // RDKit❗✔️:   const auto &oatoms = embObj.GetEmbeddedAtoms();
        // RDKit❗✔️:   for (const auto &ori : oatoms) {
        // RDKit❗✔️:     auto aid = ori.first;
        // RDKit❗✔️:     if (std::find(commAtms.begin(), commAtms.end(), aid) == commAtms.end()) {
        // RDKit❗✔️:       d_eatoms[aid] = ori.second;
        // RDKit❗✔️:       if (!ori.second.neighs.empty()) {
        // RDKit❗✔️:         if (std::find(d_attachPts.begin(), d_attachPts.end(), aid) ==
        // RDKit❗✔️:             d_attachPts.end()) {
        // RDKit❗✔️:           d_attachPts.push_back(aid);
        // RDKit❗✔️:         }
        // RDKit❗✔️:       }
        // RDKit❗✔️:     } else {
        // RDKit❗✔️:       if (ori.second.CisTransNbr >= 0) {
        // RDKit❗✔️:         d_eatoms[aid].CisTransNbr = ori.second.CisTransNbr;
        // RDKit❗✔️:         d_eatoms[aid].normal = ori.second.normal;
        // RDKit❗✔️:         d_eatoms[aid].ccw = ori.second.ccw;
        // RDKit❗✔️:       }
        // RDKit❗✔️:       if (ori.second.angle > 0.0) {
        // RDKit❗✔️:         d_eatoms[aid].angle = ori.second.angle;
        // RDKit❗✔️:         d_eatoms[aid].nbr1 = ori.second.nbr1;
        // RDKit❗✔️:         d_eatoms[aid].nbr2 = ori.second.nbr2;
        // RDKit❗✔️:       }
        // RDKit❗✔️:     }
        // RDKit❗✔️:   }
        // RDKit❗✔️:   for (auto cai : commAtms) {
        // RDKit❗✔️:     this->updateNewNeighs(cai);
        // RDKit❗✔️:   }
        // RDKit❗✔️: }
        self.check_same_topology(other)?;
        if common.is_empty() {
            return Err(FragmentError::NoCommonAtoms);
        }
        for &aid in common.iter() {
            if !self.atoms.contains_key(&aid) || !other.atoms.contains_key(&aid) {
                return Err(FragmentError::AtomNotEmbedded { atom: aid });
            }
        }
        let mut ct_case = 0;
        if common.len() == 1 {
            let aid = common[0];
            let this = &self.atoms[&aid];
            let other_atom = if this.cis_trans_nbr.is_some() {
                ct_case = 2;
                this.nbr1
            } else if other.atoms[&aid].cis_trans_nbr.is_some() {
                ct_case = 1;
                other.atoms[&aid].nbr1
            } else {
                this.nbr1
            };
            if let Some(other_aid) = other_atom {
                match ct_case {
                    1 => self.add_non_ring_atom(other_aid, aid)?,
                    _ => other.add_non_ring_atom(other_aid, aid)?,
                }
                common.push(other_aid);
            }
        }
        let transform = if common.len() == 1 {
            self.compute_one_atom_trans(common[0], other)?
        } else {
            self.compute_two_atom_trans(common[0], common[1], other)?
        };
        other.transform(transform);
        if common.len() >= 2 {
            if ct_case != 0 {
                self.reflect_if_necessary_cis_trans(other, ct_case, common[0], common[1])?;
            } else if common.len() == 2 {
                self.reflect_if_necessary_density(other, common[0], common[1])?;
            } else {
                self.reflect_if_necessary_third_pt(other, common[0], common[1], common[2])?;
            }
        }
        for (&aid, atom) in &other.atoms {
            if !common.contains(&aid) {
                // RDKit❗✔️: EmbeddedAtom &operator=(const EmbeddedAtom &other) {
                // RDKit❗✔️:   if (this == &other) {
                // RDKit❗✔️:     return *this;
                // RDKit❗✔️:   }
                // RDKit❗✔️:   loc = other.loc;
                // RDKit❗✔️:   angle = other.angle;
                // RDKit❗✔️:   nbr1 = other.nbr1;
                // RDKit❗✔️:   nbr2 = other.nbr2;
                // RDKit❗✔️:   CisTransNbr = other.CisTransNbr;
                // RDKit❗✔️:   rotDir = other.rotDir;
                // RDKit❗✔️:   normal = other.normal;
                // RDKit❗✔️:   ccw = other.ccw;
                // RDKit❗✔️:   neighs = other.neighs;
                // RDKit❗✔️:   d_density = other.d_density;
                // RDKit❗✔️:   df_fixed = other.df_fixed;
                // RDKit❗✔️:   return *this;
                // RDKit❗✔️: }
                // `operator[]` default-constructs the destination and this
                // source assignment notably does not copy `aid`.
                let mut copied = atom.clone();
                copied.aid = 0;
                self.atoms.insert(aid, copied);
                if !atom.neighs.is_empty() && !self.attachment_points.contains(&aid) {
                    self.attachment_points.push(aid);
                }
            } else {
                let target = self.atoms.get_mut(&aid).expect("common checked");
                if atom.cis_trans_nbr.is_some() {
                    target.cis_trans_nbr = atom.cis_trans_nbr;
                    target.normal = atom.normal;
                    target.ccw = atom.ccw;
                }
                if atom.angle > 0.0 {
                    target.angle = atom.angle;
                    target.nbr1 = atom.nbr1;
                    target.nbr2 = atom.nbr2;
                }
            }
        }
        for &aid in common.iter() {
            self.update_new_neighbors(aid)?;
        }
        Ok(())
    }

    fn check_same_topology(&self, other: &Self) -> Result<(), FragmentError> {
        if std::ptr::eq(self.topology, other.topology) {
            Ok(())
        } else {
            Err(FragmentError::MismatchedTopology)
        }
    }

    fn compute_one_atom_trans(
        &self,
        aid: usize,
        other: &Self,
    ) -> Result<Transform2D, FragmentError> {
        // RDKit❗✔️: RDGeom::Transform2D EmbeddedFrag::computeOneAtomTrans(
        // RDKit❗✔️:     unsigned int commAid, const EmbeddedFrag &other) {
        // RDKit❗✔️:   auto rcr = d_eatoms[commAid].loc;
        // RDKit❗✔️:   const auto &oeatm = other.GetEmbeddedAtom(commAid);
        // RDKit❗✔️:   auto ccr = oeatm.loc;
        // RDKit❗✔️:   auto onb1 = oeatm.nbr1;
        // RDKit❗✔️:   auto onb2 = oeatm.nbr2;
        // RDKit❗✔️:   CHECK_INVARIANT((onb1 >= 0) && (onb2 >= 0), "");
        // RDKit❗✔️:   auto midPt = other.GetEmbeddedAtom(onb1).loc;
        // RDKit❗✔️:   midPt += other.GetEmbeddedAtom(onb2).loc;
        // RDKit❗✔️:   midPt *= 0.5;
        // RDKit❗✔️:   auto nb1 = d_eatoms[commAid].nbr1;
        // RDKit❗✔️:   auto nb2 = d_eatoms[commAid].nbr2;
        // RDKit❗✔️:   auto nbp1 = d_eatoms[nb1].loc;
        // RDKit❗✔️:   auto nbp2 = d_eatoms[nb2].loc;
        // RDKit❗✔️:   auto ang = d_eatoms[commAid].angle;
        // RDKit❗✔️:   auto largestAngle = 2 * M_PI - ang;
        // RDKit❗✔️:   auto bpt = computeBisectPoint(rcr, largestAngle, nbp1, nbp2);
        // RDKit❗✔️:   RDGeom::Transform2D trans;
        // RDKit❗✔️:   trans.SetTransform(rcr, bpt, ccr, midPt);
        // RDKit❗✔️:   return trans;
        // RDKit❗✔️: }
        let current = self
            .atoms
            .get(&aid)
            .ok_or(FragmentError::AtomNotEmbedded { atom: aid })?;
        let incoming = other
            .atoms
            .get(&aid)
            .ok_or(FragmentError::AtomNotEmbedded { atom: aid })?;
        let i1 = incoming
            .nbr1
            .ok_or(FragmentError::NotEnoughEmbeddedNeighbors {
                atom: aid,
                count: 0,
            })?;
        let i2 = incoming
            .nbr2
            .ok_or(FragmentError::NotEnoughEmbeddedNeighbors {
                atom: aid,
                count: 1,
            })?;
        let a = other
            .atoms
            .get(&i1)
            .ok_or(FragmentError::AtomNotEmbedded { atom: i1 })?
            .loc;
        let b = other
            .atoms
            .get(&i2)
            .ok_or(FragmentError::AtomNotEmbedded { atom: i2 })?
            .loc;
        let mid = [(a[0] + b[0]) * 0.5, (a[1] + b[1]) * 0.5];
        let n1 = current
            .nbr1
            .ok_or(FragmentError::NotEnoughEmbeddedNeighbors {
                atom: aid,
                count: 0,
            })?;
        let n2 = current
            .nbr2
            .ok_or(FragmentError::NotEnoughEmbeddedNeighbors {
                atom: aid,
                count: 1,
            })?;
        let a = self
            .atoms
            .get(&n1)
            .ok_or(FragmentError::AtomNotEmbedded { atom: n1 })?
            .loc;
        let b = self
            .atoms
            .get(&n2)
            .ok_or(FragmentError::AtomNotEmbedded { atom: n2 })?
            .loc;
        let bisect = compute_bisect_point(current.loc, 2.0 * PI - current.angle, a, b);
        Ok(Transform2D::from_point_pairs(
            current.loc,
            bisect,
            incoming.loc,
            mid,
        ))
    }

    fn compute_two_atom_trans(
        &self,
        aid1: usize,
        aid2: usize,
        other: &Self,
    ) -> Result<Transform2D, FragmentError> {
        // RDKit❗✔️:     auto cid1 = commAtms[0];
        // RDKit❗✔️:     auto cid2 = commAtms[1];
        // RDKit❗✔️:     const auto &ref1 = d_eatoms.at(cid1).loc;
        // RDKit❗✔️:     const auto &ref2 = d_eatoms.at(cid2).loc;
        // RDKit❗✔️:     const auto &oth1 = embObj.GetEmbeddedAtom(cid1).loc;
        // RDKit❗✔️:     const auto &oth2 = embObj.GetEmbeddedAtom(cid2).loc;
        // RDKit❗✔️:     rtrans.SetTransform(ref1, ref2, oth1, oth2);
        let a = self
            .atoms
            .get(&aid1)
            .ok_or(FragmentError::AtomNotEmbedded { atom: aid1 })?
            .loc;
        let b = self
            .atoms
            .get(&aid2)
            .ok_or(FragmentError::AtomNotEmbedded { atom: aid2 })?
            .loc;
        let c = other
            .atoms
            .get(&aid1)
            .ok_or(FragmentError::AtomNotEmbedded { atom: aid1 })?
            .loc;
        let d = other
            .atoms
            .get(&aid2)
            .ok_or(FragmentError::AtomNotEmbedded { atom: aid2 })?
            .loc;
        Ok(Transform2D::from_point_pairs(a, b, c, d))
    }

    fn transform(&mut self, trans: Transform2D) {
        // RDKit❗✔️: void EmbeddedFrag::Transform(const RDGeom::Transform2D &trans) {
        // RDKit❗✔️:   for (auto &eri : d_eatoms) {
        // RDKit❗✔️:     eri.second.Transform(trans);
        // RDKit❗✔️:   }
        // RDKit❗✔️: }
        for atom in self.atoms.values_mut() {
            atom.transform(trans);
        }
    }

    fn compute_box(&self) -> [f64; 4] {
        // RDKit❗✔️: void EmbeddedFrag::computeBox() {
        // RDKit❗✔️:   d_px = -1.0e8;
        // RDKit❗✔️:   d_nx = 1.0e8;
        // RDKit❗✔️:   d_py = -1.0e8;
        // RDKit❗✔️:   d_ny = 1.0e8;
        // RDKit❗✔️:
        // RDKit❗✔️:   for (const auto &eri : d_eatoms) {
        // RDKit❗✔️:     const auto &loc = eri.second.loc;
        // RDKit❗✔️:     d_px = std::max(d_px, loc.x);
        // RDKit❗✔️:     d_nx = std::min(d_nx, loc.x);
        // RDKit❗✔️:     d_py = std::max(d_py, loc.y);
        // RDKit❗✔️:     d_ny = std::min(d_ny, loc.y);
        // RDKit❗✔️:   }
        // RDKit❗✔️:   d_nx *= -1.0;
        // RDKit❗✔️:   d_ny *= -1.0;
        // RDKit❗✔️: }
        // Behavior: the return order is source px/nx/py/ny, including the
        // sentinel result for an empty fragment and the final negation of
        // both minima.
        // Complexity: one ordered map traversal with constant local state,
        // matching the source without persisting otherwise unused box fields.
        let mut px: f64 = -1.0e8;
        let mut nx: f64 = 1.0e8;
        let mut py: f64 = -1.0e8;
        let mut ny: f64 = 1.0e8;
        for atom in self.atoms.values() {
            px = px.max(atom.loc[0]);
            nx = nx.min(atom.loc[0]);
            py = py.max(atom.loc[1]);
            ny = ny.min(atom.loc[1]);
        }
        [px, -nx, py, -ny]
    }

    fn canonicalize_orientation(&mut self) {
        // RDKit❗✔️: void EmbeddedFrag::canonicalizeOrientation() {
        // RDKit❗✔️:   // fix for issue 198
        // RDKit❗✔️:   // no need to canonicalize if we are dealing with a single atm
        // RDKit❗✔️:   if (d_eatoms.size() <= 1) {
        // RDKit❗✔️:     return;
        // RDKit❗✔️:   }
        // RDKit❗✔️:
        // RDKit❗✔️:   RDGeom::Point2D cent(0.0, 0.0);
        // RDKit❗✔️:   for (const auto &elem : d_eatoms) {
        // RDKit❗✔️:     cent += elem.second.loc;
        // RDKit❗✔️:   }
        // RDKit❗✔️:   cent *= (1.0 / d_eatoms.size());
        // RDKit❗✔️:
        // RDKit❗✔️:   double xx = 0.0;
        // RDKit❗✔️:   double xy = 0.0;
        // RDKit❗✔️:   double yy = 0.0;
        // RDKit❗✔️:
        // RDKit❗✔️:   // shift the center of the fragment to the origin and compute the covariance
        // RDKit❗✔️:   // matrix
        // RDKit❗✔️:   for (auto &elem : d_eatoms) {
        // RDKit❗✔️:     elem.second.loc -= cent;
        // RDKit❗✔️:     xx += (elem.second.loc.x) * (elem.second.loc.x);
        // RDKit❗✔️:     xy += (elem.second.loc.x) * (elem.second.loc.y);
        // RDKit❗✔️:     yy += (elem.second.loc.y) * (elem.second.loc.y);
        // RDKit❗✔️:   }
        // RDKit❗✔️:
        // RDKit❗✔️:   RDGeom::Point2D eig1, eig2;
        // RDKit❗✔️:   // the eigen vectors are given by
        // RDKit❗✔️:   //   (2*xy, (yy - xx) + d) and (2*xy, (yy - xx) - d)
        // RDKit❗✔️:   // where d = sqrt((xx - yy)^2 + 4*xy^2)
        // RDKit❗✔️:   auto d = (xx - yy) * (xx - yy) + 4 * xy * xy;
        // RDKit❗✔️:   d = sqrt(d);
        // RDKit❗✔️:   eig1.x = 2 * xy;
        // RDKit❗✔️:   eig1.y = (yy - xx) + d;
        // RDKit❗✔️:   if (eig1.length() <= 1e-4) {
        // RDKit❗✔️:     return;
        // RDKit❗✔️:   }
        // RDKit❗✔️:   auto eVal1 = (xx + yy + d) / 2;
        // RDKit❗✔️:   eig1.normalize();
        // RDKit❗✔️:
        // RDKit❗✔️:   eig2.x = 2 * xy;
        // RDKit❗✔️:   eig2.y = (yy - xx) - d;
        // RDKit❗✔️:   auto eVal2 = (xx + yy - d) / 2;
        // RDKit❗✔️:
        // RDKit❗✔️:   if (eig2.length() > 1e-4) {
        // RDKit❗✔️:     eig2.normalize();
        // RDKit❗✔️:
        // RDKit❗✔️:     // make sure eig1 corresponds to the larger eigenvalue:
        // RDKit❗✔️:     if (eVal2 > eVal1) {
        // RDKit❗✔️:       std::swap(eig1, eig2);
        // RDKit❗✔️:     }
        // RDKit❗✔️:   }
        // RDKit❗✔️:   // now rotate eig1 onto the X axis:
        // RDKit❗✔️:   RDGeom::Transform2D trans;
        // RDKit❗✔️:   trans.setVal(0, 0, eig1.x);
        // RDKit❗✔️:   trans.setVal(1, 0, -eig1.y);
        // RDKit❗✔️:   trans.setVal(0, 1, eig1.y);
        // RDKit❗✔️:   trans.setVal(1, 1, eig1.x);
        // RDKit❗✔️:   this->Transform(trans);
        // RDKit❗✔️: }
        // Behavior: map-order accumulation, early return after centering, the
        // 1e-4 eigenvector thresholds, eigenvalue swap and row-major rotation
        // all follow the pinned source exactly.
        // Complexity: two linear atom traversals and constant-size eigensystem
        // arithmetic, matching the source without allocation.
        if self.atoms.len() <= 1 {
            return;
        }
        let mut centroid = [0.0, 0.0];
        for atom in self.atoms.values() {
            centroid[0] += atom.loc[0];
            centroid[1] += atom.loc[1];
        }
        let scale = 1.0 / self.atoms.len() as f64;
        centroid[0] *= scale;
        centroid[1] *= scale;

        let mut xx = 0.0;
        let mut xy = 0.0;
        let mut yy = 0.0;
        for atom in self.atoms.values_mut() {
            atom.loc[0] -= centroid[0];
            atom.loc[1] -= centroid[1];
            xx += atom.loc[0] * atom.loc[0];
            xy += atom.loc[0] * atom.loc[1];
            yy += atom.loc[1] * atom.loc[1];
        }

        let d = ((xx - yy) * (xx - yy) + 4.0 * xy * xy).sqrt();
        let mut eig1 = [2.0 * xy, (yy - xx) + d];
        let eig1_length = (eig1[0] * eig1[0] + eig1[1] * eig1[1]).sqrt();
        if eig1_length <= 1.0e-4 {
            return;
        }
        let eval1 = (xx + yy + d) / 2.0;
        eig1[0] /= eig1_length;
        eig1[1] /= eig1_length;

        let mut eig2 = [2.0 * xy, (yy - xx) - d];
        let eval2 = (xx + yy - d) / 2.0;
        let eig2_length = (eig2[0] * eig2[0] + eig2[1] * eig2[1]).sqrt();
        if eig2_length > 1.0e-4 {
            eig2[0] /= eig2_length;
            eig2[1] /= eig2_length;
            if eval2 > eval1 {
                std::mem::swap(&mut eig1, &mut eig2);
            }
        }
        self.transform(Transform2D::from_linear_rows(
            [eig1[0], eig1[1]],
            [-eig1[1], eig1[0]],
        ));
    }

    fn translate(&mut self, shift: Point2) {
        // RDKit❗✔️: void Translate(const RDGeom::Point2D &shift) {
        // RDKit❗✔️:   INT_EATOM_MAP_I eari;
        // RDKit❗✔️:   for (eari = d_eatoms.begin(); eari != d_eatoms.end(); eari++) {
        // RDKit❗✔️:     eari->second.loc += shift;
        // RDKit❗✔️:   }
        // RDKit❗✔️: }
        // Behavior: translation changes locations only; normals and all
        // non-coordinate embedded-atom state remain untouched.
        // Complexity: one ordered linear traversal with no allocation.
        for atom in self.atoms.values_mut() {
            atom.loc[0] += shift[0];
            atom.loc[1] += shift[1];
        }
    }

    fn reflect(&mut self, a: Point2, b: Point2) {
        // RDKit❗✔️: void EmbeddedFrag::Reflect(const RDGeom::Point2D &loc1,
        // RDKit❗✔️:                            const RDGeom::Point2D &loc2) {
        // RDKit❗✔️:   for (auto &ei : d_eatoms) {
        // RDKit❗✔️:     ei.second.Reflect(loc1, loc2);
        // RDKit❗✔️:   }
        // RDKit❗✔️: }
        for atom in self.atoms.values_mut() {
            atom.reflect(a, b);
        }
    }

    fn reflect_if_necessary_third_pt(
        &self,
        other: &mut Self,
        aid1: usize,
        aid2: usize,
        aid3: usize,
    ) -> Result<(), FragmentError> {
        // RDKit❗✔️: void EmbeddedFrag::reflectIfNecessaryThirdPt(EmbeddedFrag &embFrag,
        // RDKit❗✔️:                                              unsigned int aid1,
        // RDKit❗✔️:                                              unsigned int aid2,
        // RDKit❗✔️:                                              unsigned int aid3) {
        // RDKit❗✔️:   const auto &pt1 = d_eatoms[aid1].loc;
        // RDKit❗✔️:   const auto &pt2 = d_eatoms[aid2].loc;
        // RDKit❗✔️:   auto normal = pt2;
        // RDKit❗✔️:   normal -= pt1;
        // RDKit❗✔️:   normal.rotate90();
        // RDKit❗✔️:   const auto oth3 = embFrag.GetEmbeddedAtom(aid3).loc - pt1;
        // RDKit❗✔️:   const auto pt3 = d_eatoms[aid3].loc - pt1;
        // RDKit❗✔️:   auto dot1 = normal.dotProduct(pt3);
        // RDKit❗✔️:   auto dot2 = normal.dotProduct(oth3);
        // RDKit❗✔️:   if (dot1 * dot2 < 0.0) {
        // RDKit❗✔️:     embFrag.Reflect(pt1, pt2);
        // RDKit❗✔️:   }
        // RDKit❗✔️: }
        let p1 = self
            .atoms
            .get(&aid1)
            .ok_or(FragmentError::AtomNotEmbedded { atom: aid1 })?
            .loc;
        let p2 = self
            .atoms
            .get(&aid2)
            .ok_or(FragmentError::AtomNotEmbedded { atom: aid2 })?
            .loc;
        let p3 = self
            .atoms
            .get(&aid3)
            .ok_or(FragmentError::AtomNotEmbedded { atom: aid3 })?
            .loc;
        let q3 = other
            .atoms
            .get(&aid3)
            .ok_or(FragmentError::AtomNotEmbedded { atom: aid3 })?
            .loc;
        let normal = [-(p2[1] - p1[1]), p2[0] - p1[0]];
        let dot1 = normal[0] * (p3[0] - p1[0]) + normal[1] * (p3[1] - p1[1]);
        let dot2 = normal[0] * (q3[0] - p1[0]) + normal[1] * (q3[1] - p1[1]);
        if dot1 * dot2 < 0.0 {
            other.reflect(p1, p2);
        }
        Ok(())
    }

    fn reflect_if_necessary_density(
        &self,
        other: &mut Self,
        aid1: usize,
        aid2: usize,
    ) -> Result<(), FragmentError> {
        // RDKit❗✔️: void EmbeddedFrag::reflectIfNecessaryDensity(EmbeddedFrag &embFrag,
        // RDKit❗✔️:                                              unsigned int aid1,
        // RDKit❗✔️:                                              unsigned int aid2) {
        // RDKit❗✔️:   const auto &pin1 = d_eatoms[aid1].loc;
        // RDKit❗✔️:   const auto &pin2 = d_eatoms[aid2].loc;
        // RDKit❗✔️:   double densityNormal = 0.0;
        // RDKit❗✔️:   double densityReflect = 0.0;
        // RDKit❗✔️:   for (const auto &oci : embFrag.GetEmbeddedAtoms()) {
        // RDKit❗✔️:     if (d_eatoms.find(oci.first) == d_eatoms.end()) {
        // RDKit❗✔️:       auto loc1 = oci.second.loc;
        // RDKit❗✔️:       auto rloc1 = reflectPoint(loc1, pin1, pin2);
        // RDKit❗✔️:       for (const auto &tci : d_eatoms) {
        // RDKit❗✔️:         auto t1 = tci.second.loc;
        // RDKit❗✔️:         t1 -= loc1;
        // RDKit❗✔️:         auto td = t1.length();
        // RDKit❗✔️:         auto rt1 = tci.second.loc;
        // RDKit❗✔️:         rt1 -= rloc1;
        // RDKit❗✔️:         auto rtd = rt1.length();
        // RDKit❗✔️:         if (td > 1.0e-3) {
        // RDKit❗✔️:           densityNormal += (1.0 / td);
        // RDKit❗✔️:         } else {
        // RDKit❗✔️:           densityNormal += 1000.0;
        // RDKit❗✔️:         }
        // RDKit❗✔️:         if (rtd > 1.0e-3) {
        // RDKit❗✔️:           densityReflect += (1.0 / rtd);
        // RDKit❗✔️:         } else {
        // RDKit❗✔️:           densityReflect += 1000.0;
        // RDKit❗✔️:         }
        // RDKit❗✔️:       }
        // RDKit❗✔️:     }
        // RDKit❗✔️:   }
        // RDKit❗✔️:   if (densityNormal - densityReflect > 1.0e-4) {
        // RDKit❗✔️:     embFrag.Reflect(pin1, pin2);
        // RDKit❗✔️:   }
        // RDKit❗✔️: }
        let a = self
            .atoms
            .get(&aid1)
            .ok_or(FragmentError::AtomNotEmbedded { atom: aid1 })?
            .loc;
        let b = self
            .atoms
            .get(&aid2)
            .ok_or(FragmentError::AtomNotEmbedded { atom: aid2 })?
            .loc;
        let mut normal_density = 0.0;
        let mut reflected_density = 0.0;
        for (&aid, atom) in &other.atoms {
            if !self.atoms.contains_key(&aid) {
                let point = atom.loc;
                let reflected = reflect_point(point, a, b);
                for target in self.atoms.values() {
                    let d = (target.loc[0] - point[0]).hypot(target.loc[1] - point[1]);
                    let r = (target.loc[0] - reflected[0]).hypot(target.loc[1] - reflected[1]);
                    normal_density += if d > 1e-3 { 1.0 / d } else { 1000.0 };
                    reflected_density += if r > 1e-3 { 1.0 / r } else { 1000.0 };
                }
            }
        }
        if normal_density - reflected_density > 1e-4 {
            other.reflect(a, b);
        }
        Ok(())
    }

    fn reflect_if_necessary_cis_trans(
        &self,
        other: &mut Self,
        ct_case: u8,
        aid1: usize,
        aid2: usize,
    ) -> Result<(), FragmentError> {
        // RDKit❗✔️: void EmbeddedFrag::reflectIfNecessaryCisTrans(EmbeddedFrag &embFrag,
        // RDKit❗✔️:                                               unsigned int ctCase,
        // RDKit❗✔️:                                               unsigned int aid1,
        // RDKit❗✔️:                                               unsigned int aid2) {
        // RDKit❗✔️:   const auto &p1Loc = d_eatoms[aid1].loc;
        // RDKit❗✔️:   RDGeom::Point2D rAtmLoc, p1norm;
        // RDKit❗✔️:   if (ctCase == 1) {
        // RDKit❗✔️:     p1norm = embFrag.d_eatoms[aid1].normal;
        // RDKit❗✔️:     auto ringAtm = embFrag.d_eatoms[aid1].CisTransNbr;
        // RDKit❗✔️:     if (d_eatoms.find(ringAtm) != d_eatoms.end()) {
        // RDKit❗✔️:       rAtmLoc = d_eatoms[ringAtm].loc;
        // RDKit❗✔️:     } else {
        // RDKit❗✔️:       BOOST_LOG(rdWarningLog) << "Warning: stereochemistry around double bond "
        // RDKit❗✔️:                                  "may be incorrect in depiction."
        // RDKit❗✔️:                               << std::endl;
        // RDKit❗✔️:       return;
        // RDKit❗✔️:     }
        // RDKit❗✔️:   } else {
        // RDKit❗✔️:     p1norm = d_eatoms[aid1].normal;
        // RDKit❗✔️:     auto ringAtm = d_eatoms[aid1].CisTransNbr;
        // RDKit❗✔️:     rAtmLoc = embFrag.d_eatoms[ringAtm].loc;
        // RDKit❗✔️:   }
        // RDKit❗✔️:   rAtmLoc -= p1Loc;
        // RDKit❗✔️:   auto dot = rAtmLoc.dotProduct(p1norm);
        // RDKit❗✔️:   auto p2Loc = d_eatoms[aid2].loc;
        // RDKit❗✔️:   if (dot < 0.0) {
        // RDKit❗✔️:     embFrag.Reflect(p1Loc, p2Loc);
        // RDKit❗✔️:   }
        // RDKit❗✔️: }
        let origin = self
            .atoms
            .get(&aid1)
            .ok_or(FragmentError::AtomNotEmbedded { atom: aid1 })?
            .loc;
        let (normal, ring) = if ct_case == 1 {
            let atom = other
                .atoms
                .get(&aid1)
                .ok_or(FragmentError::AtomNotEmbedded { atom: aid1 })?;
            (atom.normal, atom.cis_trans_nbr)
        } else {
            let atom = self
                .atoms
                .get(&aid1)
                .ok_or(FragmentError::AtomNotEmbedded { atom: aid1 })?;
            (atom.normal, atom.cis_trans_nbr)
        };
        let ring = ring.ok_or(FragmentError::NotEnoughEmbeddedNeighbors {
            atom: aid1,
            count: 0,
        })?;
        let ring_loc = if ct_case == 1 {
            // The pinned warning path deliberately leaves the transformed
            // fragment unchanged when the ring neighbor is absent.
            let Some(atom) = self.atoms.get(&ring) else {
                return Ok(());
            };
            atom.loc
        } else {
            other
                .atoms
                .get(&ring)
                .ok_or(FragmentError::AtomNotEmbedded { atom: ring })?
                .loc
        };
        let dot = (ring_loc[0] - origin[0]) * normal[0] + (ring_loc[1] - origin[1]) * normal[1];
        if dot < 0.0 {
            let second = self
                .atoms
                .get(&aid2)
                .ok_or(FragmentError::AtomNotEmbedded { atom: aid2 })?
                .loc;
            other.reflect(origin, second);
        }
        Ok(())
    }
}

fn compute_sub_angle(degree: usize, hybridization: Hybridization) -> f64 {
    // RDKit❗✔️: inline double computeSubAngle(unsigned int degree,
    // RDKit❗✔️:                               RDKit::Atom::HybridizationType htype) {
    // RDKit❗✔️:   double angle = M_PI;
    // RDKit❗✔️:   switch (htype) {
    // RDKit❗✔️:     case RDKit::Atom::UNSPECIFIED:
    // RDKit❗✔️:     case RDKit::Atom::SP3:
    // RDKit❗✔️:       if (degree == 4) {
    // RDKit❗✔️:         angle = M_PI / 2;
    // RDKit❗✔️:       } else {
    // RDKit❗✔️:         angle = 2 * M_PI / 3;
    // RDKit❗✔️:       }
    // RDKit❗✔️:       break;
    // RDKit❗✔️:     case RDKit::Atom::SP2:
    // RDKit❗✔️:       angle = 2 * M_PI / 3;
    // RDKit❗✔️:       break;
    // RDKit❗✔️:     default:
    // RDKit❗✔️:       angle = 2. * M_PI / degree;
    // RDKit❗✔️:   }
    // RDKit❗✔️:   return angle;
    // RDKit❗✔️: }
    match hybridization {
        Hybridization::Unspecified | Hybridization::Sp3 if degree == 4 => PI / 2.0,
        Hybridization::Unspecified | Hybridization::Sp3 | Hybridization::Sp2 => 2.0 * PI / 3.0,
        _ => 2.0 * PI / degree as f64,
    }
}

fn normalize(point: Point2) -> Result<Point2, FragmentError> {
    // RDKit❗✔️: void normalize() override {
    // RDKit❗✔️:   double ln = this->length();
    // RDKit❗✔️:   if (ln < zero_tolerance) {
    // RDKit❗✔️:     throw std::runtime_error("Cannot normalize a zero length vector");
    // RDKit❗✔️:   }
    // RDKit❗✔️:   x /= ln;
    // RDKit❗✔️:   y /= ln;
    // RDKit❗✔️: }
    let length = (point[0] * point[0] + point[1] * point[1]).sqrt();
    // `RDGeneral/Numerics/Vector.h::zero_tolerance` is 1.e-16 in the pin.
    if length < 1e-16 {
        return Err(FragmentError::CoincidentPoints);
    }
    Ok([point[0] / length, point[1] / length])
}

fn compute_normal(center: Point2, other: Point2) -> Result<Point2, FragmentError> {
    // RDKit❗✔️: inline RDGeom::Point2D computeNormal(const RDGeom::Point2D &center,
    // RDKit❗✔️:                                      const RDGeom::Point2D &other) {
    // RDKit❗✔️:   auto res = other - center;
    // RDKit❗✔️:   res.normalize();
    // RDKit❗✔️:   std::swap(res.x, res.y);
    // RDKit❗✔️:   res.x *= -1;
    // RDKit❗✔️:   return res;
    // RDKit❗✔️: }
    let vector = normalize([other[0] - center[0], other[1] - center[1]])?;
    Ok([-vector[1], vector[0]])
}

fn compute_angle(center: Point2, loc1: Point2, loc2: Point2) -> Result<f64, FragmentError> {
    // RDKit❗✔️: inline double computeAngle(const RDGeom::Point2D &center,
    // RDKit❗✔️:                            const RDGeom::Point2D &loc1,
    // RDKit❗✔️:                            const RDGeom::Point2D &loc2) {
    // RDKit❗✔️:   auto v1 = loc1 - center;
    // RDKit❗✔️:   auto v2 = loc2 - center;
    // RDKit❗✔️:   return v1.angleTo(v2);
    // RDKit❗✔️: }
    // RDKit❗✔️: double angleTo(const Point2D &other) const {
    // RDKit❗✔️:   auto t1 = *this;
    // RDKit❗✔️:   auto t2 = other;
    // RDKit❗✔️:   t1.normalize();
    // RDKit❗✔️:   t2.normalize();
    // RDKit❗✔️:   double dotProd = t1.dotProduct(t2);
    // RDKit❗✔️:   if (dotProd < -1.0) {
    // RDKit❗✔️:     dotProd = -1.0;
    // RDKit❗✔️:   } else if (dotProd > 1.0) {
    // RDKit❗✔️:     dotProd = 1.0;
    // RDKit❗✔️:   }
    // RDKit❗✔️:   return acos(dotProd);
    // RDKit❗✔️: }
    let a = normalize([loc1[0] - center[0], loc1[1] - center[1]])?;
    let b = normalize([loc2[0] - center[0], loc2[1] - center[1]])?;
    Ok((a[0] * b[0] + a[1] * b[1]).clamp(-1.0, 1.0).acos())
}

fn rotation_dir(center: Point2, loc1: Point2, loc2: Point2, remaining_angle: f64) -> i32 {
    // RDKit❗✔️: inline int rotationDir(const RDGeom::Point2D &center,
    // RDKit❗✔️:                        const RDGeom::Point2D &loc1, const RDGeom::Point2D &loc2,
    // RDKit❗✔️:                        double remAngle) {
    // RDKit❗✔️:   auto pt1 = loc1 - center;
    // RDKit❗✔️:   auto pt2 = loc2 - center;
    // RDKit❗✔️:   auto cross = pt1.x * pt2.y - pt1.y * pt2.x;
    // RDKit❗✔️:   auto diffAngle = M_PI - remAngle;
    // RDKit❗✔️:   cross *= diffAngle;
    // RDKit❗✔️:   if (cross >= 0.0) {
    // RDKit❗✔️:     return -1;
    // RDKit❗✔️:   } else {
    // RDKit❗✔️:     return 1;
    // RDKit❗✔️:   }
    // RDKit❗✔️: }
    let a = [loc1[0] - center[0], loc1[1] - center[1]];
    let b = [loc2[0] - center[0], loc2[1] - center[1]];
    if (a[0] * b[1] - a[1] * b[0]) * (PI - remaining_angle) >= 0.0 {
        -1
    } else {
        1
    }
}
