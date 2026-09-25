//! Source-shaped seed fragments for non-tetrahedral stereocentres.

use cosmolkit_core::{
    RingInfo, non_tetrahedral_across_ligand, non_tetrahedral_ideal_angle,
    trigonal_bipyramidal_axial_ligand,
};
use cosmolkit_model::{ChiralTag, TopologyBlock};

use crate::embedded_frag::{EmbeddedFrag, FragmentError};
use crate::geometry::{BOND_LEN, PointMap, atom_depict_rank};

const ISQRT2: f64 = 0.707107;
const SQRT3_2: f64 = 0.866025;

fn ranked_neighbors(topology: &TopologyBlock, centre: usize, ranks: &[i32]) -> Vec<usize> {
    // RDKit❗✔️: std::vector<const RDKit::Atom *> getRankedAtomNeighbors(
    // RDKit❗✔️:     const RDKit::ROMol &mol, const RDKit::Atom *atom,
    // RDKit❗✔️:     const std::vector<int> &atomRanks) {
    // RDKit❗✔️:   std::vector<const RDKit::Atom *> nbrs;
    // RDKit❗✔️:   for (auto nbr : mol.atomNeighbors(atom)) {
    // RDKit❗✔️:     nbrs.push_back(nbr);
    // RDKit❗✔️:   }
    // RDKit❗✔️:   std::sort(nbrs.begin(), nbrs.end(),
    // RDKit❗✔️:             [&atomRanks](const auto e1, const auto e2) {
    // RDKit❗✔️:               return atomRanks[e1->getIdx()] < atomRanks[e2->getIdx()];
    // RDKit❗✔️:             });
    // RDKit❗✔️:   return nbrs;
    // RDKit❗✔️: }
    // Behavior: ties retain adjacency insertion order in this model; source
    // std::sort does not specify their order, so exact tie parity is unproved.
    // Complexity: one degree-sized sort and allocation, matching the source.
    let mut neighbors: Vec<_> = topology
        .adjacency
        .neighbors_of(centre)
        .iter()
        .map(|neighbor| neighbor.atom_index)
        .collect();
    neighbors.sort_by_key(|&neighbor| ranks[neighbor]);
    neighbors
}

fn square_planar_points(
    topology: &TopologyBlock,
    centre: usize,
    ranks: &[i32],
) -> Result<PointMap, FragmentError> {
    // RDKit❗✔️: void embedSquarePlanar(const RDKit::ROMol &mol, const RDKit::Atom *atom,
    // RDKit❗✔️:                        std::list<EmbeddedFrag> &efrags,
    // RDKit❗✔️:                        const std::vector<int> &atomRanks) {
    // RDKit❗✔️:   static const RDGeom::Point2D idealPoints[] = {
    // RDKit❗✔️:       RDGeom::Point2D(ISQRT2 * BOND_LEN, ISQRT2 * BOND_LEN),
    // RDKit❗✔️:       RDGeom::Point2D(ISQRT2 * BOND_LEN, -ISQRT2 * BOND_LEN),
    // RDKit❗✔️:       RDGeom::Point2D(-ISQRT2 * BOND_LEN, -ISQRT2 * BOND_LEN),
    // RDKit❗✔️:       RDGeom::Point2D(-ISQRT2 * BOND_LEN, ISQRT2 * BOND_LEN),
    // RDKit❗✔️:   };
    // RDKit❗✔️:   PRECONDITION(atom, "bad atom");
    // RDKit❗✔️:   if (atom->getChiralTag() != RDKit::Atom::ChiralType::CHI_SQUAREPLANAR) {
    // RDKit❗✔️:     return;
    // RDKit❗✔️:   }
    // RDKit❗✔️:   auto nbrs = getRankedAtomNeighbors(mol, atom, atomRanks);
    // RDKit❗✔️:   RDGeom::INT_POINT2D_MAP coordMap;
    // RDKit❗✔️:   coordMap[atom->getIdx()] = RDGeom::Point2D(0., 0.);
    // RDKit❗✔️:   coordMap[nbrs[0]->getIdx()] = idealPoints[0];
    // RDKit❗✔️:   bool q2Full = false;
    // RDKit❗✔️:   for (const auto nbr : nbrs) {
    // RDKit❗✔️:     if (nbr == nbrs.front()) {
    // RDKit❗✔️:       continue;
    // RDKit❗✔️:     }
    // RDKit❗✔️:     auto angle =
    // RDKit❗✔️:         RDKit::Chirality::getIdealAngleBetweenLigands(atom, nbrs.front(), nbr);
    // RDKit❗✔️:     if (fabs(angle - 180) < 0.1) {
    // RDKit❗✔️:       coordMap[nbr->getIdx()] = idealPoints[2];
    // RDKit❗✔️:     } else {
    // RDKit❗✔️:       if (!q2Full) {
    // RDKit❗✔️:         coordMap[nbr->getIdx()] = idealPoints[1];
    // RDKit❗✔️:         q2Full = true;
    // RDKit❗✔️:       } else {
    // RDKit❗✔️:         coordMap[nbr->getIdx()] = idealPoints[3];
    // RDKit❗✔️:       }
    // RDKit❗✔️:     }
    // RDKit❗✔️:   }
    // RDKit❗✔️:   efrags.emplace_back(&mol, coordMap);
    // RDKit❗✔️: }
    // Behavior: no-neighbor source indexing is undefined; return a typed error.
    // Complexity: one sort and bounded coordinate insertion, as in source.
    let ideal = [
        [ISQRT2 * BOND_LEN, ISQRT2 * BOND_LEN],
        [ISQRT2 * BOND_LEN, -ISQRT2 * BOND_LEN],
        [-ISQRT2 * BOND_LEN, -ISQRT2 * BOND_LEN],
        [-ISQRT2 * BOND_LEN, ISQRT2 * BOND_LEN],
    ];
    let neighbors = ranked_neighbors(topology, centre, ranks);
    let first = *neighbors
        .first()
        .ok_or(FragmentError::NonTetrahedralNoLigand { centre })?;
    let mut points = PointMap::from([(centre, [0.0, 0.0]), (first, ideal[0])]);
    let mut q2_full = false;
    for &neighbor in neighbors.iter().skip(1) {
        let angle = non_tetrahedral_ideal_angle(topology, centre, first, neighbor);
        let index = if (angle - 180.0).abs() < 0.1 {
            2
        } else if !q2_full {
            q2_full = true;
            1
        } else {
            3
        };
        points.insert(neighbor, ideal[index]);
    }
    Ok(points)
}

fn trigonal_bipyramidal_points(
    topology: &TopologyBlock,
    centre: usize,
    ranks: &[i32],
) -> Result<PointMap, FragmentError> {
    // RDKit❗✔️: void embedTBP(const RDKit::ROMol &mol, const RDKit::Atom *atom,
    // RDKit❗✔️:               std::list<EmbeddedFrag> &efrags,
    // RDKit❗✔️:               const std::vector<int> &atomRanks) {
    // RDKit❗✔️:   static const RDGeom::Point2D idealPoints[] = {
    // RDKit❗✔️:       RDGeom::Point2D(0, BOND_LEN),                        // axial
    // RDKit❗✔️:       RDGeom::Point2D(0, -BOND_LEN),                       // axial
    // RDKit❗✔️:       RDGeom::Point2D(-SQRT3_2 * BOND_LEN, BOND_LEN / 2),  // equatorial
    // RDKit❗✔️:       RDGeom::Point2D(-SQRT3_2 * BOND_LEN,
    // RDKit❗✔️:                       -BOND_LEN / 2),  // equatorial
    // RDKit❗✔️:       RDGeom::Point2D(BOND_LEN, 0),    // equatorial
    // RDKit❗✔️:   };
    // RDKit❗✔️:   PRECONDITION(atom, "bad atom");
    // RDKit❗✔️:   if (atom->getChiralTag() !=
    // RDKit❗✔️:       RDKit::Atom::ChiralType::CHI_TRIGONALBIPYRAMIDAL) {
    // RDKit❗✔️:     return;
    // RDKit❗✔️:   }
    // RDKit❗✔️:   auto nbrs = getRankedAtomNeighbors(mol, atom, atomRanks);
    // RDKit❗✔️:   RDGeom::INT_POINT2D_MAP coordMap;
    // RDKit❗✔️:   coordMap[atom->getIdx()] = RDGeom::Point2D(0., 0.);
    // RDKit❗✔️:   const RDKit::Atom *axial1 =
    // RDKit❗✔️:       RDKit::Chirality::getTrigonalBipyramidalAxialAtom(atom);
    // RDKit❗✔️:   const RDKit::Atom *axial2 =
    // RDKit❗✔️:       RDKit::Chirality::getTrigonalBipyramidalAxialAtom(atom, -1);
    // RDKit❗✔️:   if (axial1) {
    // RDKit❗✔️:     coordMap[axial1->getIdx()] = idealPoints[0];
    // RDKit❗✔️:   }
    // RDKit❗✔️:   if (axial2) {
    // RDKit❗✔️:     coordMap[axial2->getIdx()] = idealPoints[1];
    // RDKit❗✔️:   }
    // RDKit❗✔️:   unsigned whichEq = 2;
    // RDKit❗✔️:   for (const auto nbr : nbrs) {
    // RDKit❗✔️:     if (nbr != axial1 && nbr != axial2) {
    // RDKit❗✔️:       coordMap[nbr->getIdx()] = idealPoints[whichEq++];
    // RDKit❗✔️:     }
    // RDKit❗✔️:   }
    // RDKit❗✔️:   efrags.emplace_back(&mol, coordMap);
    // RDKit❗✔️: }
    // Behavior: excess equatorial indexing is source-undefined, typed error here.
    // Complexity: one sort and bounded coordinate insertion, as in source.
    let ideal = [
        [0.0, BOND_LEN],
        [0.0, -BOND_LEN],
        [-SQRT3_2 * BOND_LEN, BOND_LEN / 2.0],
        [-SQRT3_2 * BOND_LEN, -BOND_LEN / 2.0],
        [BOND_LEN, 0.0],
    ];
    let axial1 = trigonal_bipyramidal_axial_ligand(topology, centre, 1);
    let axial2 = trigonal_bipyramidal_axial_ligand(topology, centre, -1);
    let mut points = PointMap::from([(centre, [0.0, 0.0])]);
    if let Some(index) = axial1 {
        points.insert(index, ideal[0]);
    }
    if let Some(index) = axial2 {
        points.insert(index, ideal[1]);
    }
    let mut equatorial = 2;
    for neighbor in ranked_neighbors(topology, centre, ranks) {
        if Some(neighbor) != axial1 && Some(neighbor) != axial2 {
            let point = ideal
                .get(equatorial)
                .ok_or(FragmentError::NonTetrahedralLigandOverflow { centre })?;
            points.insert(neighbor, *point);
            equatorial += 1;
        }
    }
    Ok(points)
}

fn octahedral_points(topology: &TopologyBlock, centre: usize, ranks: &[i32]) -> PointMap {
    // RDKit❗✔️: void embedOctahedral(const RDKit::ROMol &mol, const RDKit::Atom *atom,
    // RDKit❗✔️:                      std::list<EmbeddedFrag> &efrags,
    // RDKit❗✔️:                      const std::vector<int> &atomRanks) {
    // RDKit❗✔️:   static const RDGeom::Point2D idealPoints[] = {
    // RDKit❗✔️:       RDGeom::Point2D(0, BOND_LEN),                         // axial
    // RDKit❗✔️:       RDGeom::Point2D(0, -BOND_LEN),                        // axial
    // RDKit❗✔️:       RDGeom::Point2D(SQRT3_2 * BOND_LEN, BOND_LEN / 2),    // equatorial
    // RDKit❗✔️:       RDGeom::Point2D(SQRT3_2 * BOND_LEN, -BOND_LEN / 2),   // equatorial
    // RDKit❗✔️:       RDGeom::Point2D(-SQRT3_2 * BOND_LEN, -BOND_LEN / 2),  // equatorial
    // RDKit❗✔️:       RDGeom::Point2D(-SQRT3_2 * BOND_LEN, BOND_LEN / 2),   // equatorial
    // RDKit❗✔️:   };
    // RDKit❗✔️:   PRECONDITION(atom, "bad atom");
    // RDKit❗✔️:   if (atom->getChiralTag() != RDKit::Atom::ChiralType::CHI_OCTAHEDRAL) {
    // RDKit❗✔️:     return;
    // RDKit❗✔️:   }
    // RDKit❗✔️:   auto nbrs = getRankedAtomNeighbors(mol, atom, atomRanks);
    // RDKit❗✔️:   RDGeom::INT_POINT2D_MAP coordMap;
    // RDKit❗✔️:   coordMap[atom->getIdx()] = RDGeom::Point2D(0., 0.);
    // RDKit❗✔️:   const RDKit::Atom *axial1 = nullptr;
    // RDKit❗✔️:   const RDKit::Atom *axial2 = nullptr;
    // RDKit❗✔️:   for (auto i = 0u; i < nbrs.size(); ++i) {
    // RDKit❗✔️:     bool all90 = true;
    // RDKit❗✔️:     for (auto j = i + 1; j < nbrs.size(); ++j) {
    // RDKit❗✔️:       if (fabs(RDKit::Chirality::getIdealAngleBetweenLigands(atom, nbrs[i],
    // RDKit❗✔️:                                                              nbrs[j]) -
    // RDKit❗✔️:                180) < 0.1) {
    // RDKit❗✔️:         axial1 = nbrs[i];
    // RDKit❗✔️:         axial2 = nbrs[j];
    // RDKit❗✔️:         all90 = false;
    // RDKit❗✔️:         break;
    // RDKit❗✔️:       } else if (fabs(RDKit::Chirality::getIdealAngleBetweenLigands(
    // RDKit❗✔️:                           atom, nbrs[i], nbrs[j]) -
    // RDKit❗✔️:                       90) > 0.1) {
    // RDKit❗✔️:         all90 = false;
    // RDKit❗✔️:       }
    // RDKit❗✔️:     }
    // RDKit❗✔️:     if (all90) {
    // RDKit❗✔️:       axial1 = nbrs[i];
    // RDKit❗✔️:     }
    // RDKit❗✔️:     if (axial1) {
    // RDKit❗✔️:       break;
    // RDKit❗✔️:     }
    // RDKit❗✔️:   }
    // RDKit❗✔️:   if (axial1) {
    // RDKit❗✔️:     coordMap[axial1->getIdx()] = idealPoints[0];
    // RDKit❗✔️:   }
    // RDKit❗✔️:   if (axial2) {
    // RDKit❗✔️:     coordMap[axial2->getIdx()] = idealPoints[1];
    // RDKit❗✔️:   }
    // RDKit❗✔️:   const RDKit::Atom *refEqAtom1 = nullptr;
    // RDKit❗✔️:   const RDKit::Atom *refEqAtom2 = nullptr;
    // RDKit❗✔️:   for (const auto nbr : nbrs) {
    // RDKit❗✔️:     if (nbr != axial1 && nbr != axial2) {
    // RDKit❗✔️:       if (!refEqAtom1) {
    // RDKit❗✔️:         refEqAtom1 = nbr;
    // RDKit❗✔️:         coordMap[nbr->getIdx()] = idealPoints[2];
    // RDKit❗✔️:         refEqAtom2 = RDKit::Chirality::getChiralAcrossAtom(atom, nbr);
    // RDKit❗✔️:         if (refEqAtom2) {
    // RDKit❗✔️:           coordMap[refEqAtom2->getIdx()] = idealPoints[4];
    // RDKit❗✔️:         }
    // RDKit❗✔️:       } else {
    // RDKit❗✔️:         if (nbr == refEqAtom2 || nbr == refEqAtom1) {
    // RDKit❗✔️:           continue;
    // RDKit❗✔️:         }
    // RDKit❗✔️:         coordMap[nbr->getIdx()] = idealPoints[3];
    // RDKit❗✔️:         const auto acrossAtom2 =
    // RDKit❗✔️:             RDKit::Chirality::getChiralAcrossAtom(atom, nbr);
    // RDKit❗✔️:         if (acrossAtom2) {
    // RDKit❗✔️:           coordMap[acrossAtom2->getIdx()] = idealPoints[5];
    // RDKit❗✔️:         }
    // RDKit❗✔️:         break;
    // RDKit❗✔️:       }
    // RDKit❗✔️:     }
    // RDKit❗✔️:   }
    // RDKit❗✔️:   efrags.emplace_back(&mol, coordMap);
    // RDKit❗✔️: }
    // Behavior: incomplete ligand sets leave unplaced neighbors for later growth.
    // Complexity: bounded nested neighbor scans and map inserts, as source.
    let ideal = [
        [0.0, BOND_LEN],
        [0.0, -BOND_LEN],
        [SQRT3_2 * BOND_LEN, BOND_LEN / 2.0],
        [SQRT3_2 * BOND_LEN, -BOND_LEN / 2.0],
        [-SQRT3_2 * BOND_LEN, -BOND_LEN / 2.0],
        [-SQRT3_2 * BOND_LEN, BOND_LEN / 2.0],
    ];
    let neighbors = ranked_neighbors(topology, centre, ranks);
    let mut points = PointMap::from([(centre, [0.0, 0.0])]);
    let (mut axial1, mut axial2) = (None, None);
    for (i, &first) in neighbors.iter().enumerate() {
        let mut all_90 = true;
        for &second in neighbors.iter().skip(i + 1) {
            let angle = non_tetrahedral_ideal_angle(topology, centre, first, second);
            if (angle - 180.0).abs() < 0.1 {
                axial1 = Some(first);
                axial2 = Some(second);
                all_90 = false;
                break;
            } else if (angle - 90.0).abs() > 0.1 {
                all_90 = false;
            }
        }
        if all_90 {
            axial1 = Some(first);
        }
        if axial1.is_some() {
            break;
        }
    }
    if let Some(index) = axial1 {
        points.insert(index, ideal[0]);
    }
    if let Some(index) = axial2 {
        points.insert(index, ideal[1]);
    }
    let (mut first_equatorial, mut second_equatorial) = (None, None);
    for neighbor in neighbors {
        if Some(neighbor) == axial1 || Some(neighbor) == axial2 {
            continue;
        }
        if first_equatorial.is_none() {
            first_equatorial = Some(neighbor);
            points.insert(neighbor, ideal[2]);
            second_equatorial = non_tetrahedral_across_ligand(topology, centre, neighbor);
            if let Some(index) = second_equatorial {
                points.insert(index, ideal[4]);
            }
        } else {
            if Some(neighbor) == second_equatorial || Some(neighbor) == first_equatorial {
                continue;
            }
            points.insert(neighbor, ideal[3]);
            if let Some(index) = non_tetrahedral_across_ligand(topology, centre, neighbor) {
                points.insert(index, ideal[5]);
            }
            break;
        }
    }
    points
}

pub(crate) fn embed_nontetrahedral_stereo<'a>(
    topology: &'a TopologyBlock,
    rings: &'a RingInfo,
) -> Result<Vec<EmbeddedFrag<'a>>, FragmentError> {
    // RDKit❗✔️: void embedNontetrahedralStereo(const RDKit::ROMol &mol,
    // RDKit❗✔️:                                std::list<EmbeddedFrag> &efrags,
    // RDKit❗✔️:                                const std::vector<int> &atomRanks) {
    // RDKit❗✔️:   boost::dynamic_bitset<> consider(mol.getNumAtoms());
    // RDKit❗✔️:   for (const auto atm : mol.atoms()) {
    // RDKit❗✔️:     if (RDKit::Chirality::hasNonTetrahedralStereo(atm)) {
    // RDKit❗✔️:       consider[atm->getIdx()] = 1;
    // RDKit❗✔️:     }
    // RDKit❗✔️:   }
    // RDKit❗✔️:   if (consider.empty()) {
    // RDKit❗✔️:     return;
    // RDKit❗✔️:   }
    // RDKit❗✔️:   for (const auto atm : mol.atoms()) {
    // RDKit❗✔️:     if (!consider[atm->getIdx()]) {
    // RDKit❗✔️:       continue;
    // RDKit❗✔️:     }
    // RDKit❗✔️:     switch (atm->getChiralTag()) {
    // RDKit❗✔️:       case RDKit::Atom::ChiralType::CHI_SQUAREPLANAR:
    // RDKit❗✔️:         embedSquarePlanar(mol, atm, efrags, atomRanks);
    // RDKit❗✔️:         break;
    // RDKit❗✔️:       case RDKit::Atom::ChiralType::CHI_TRIGONALBIPYRAMIDAL:
    // RDKit❗✔️:         embedTBP(mol, atm, efrags, atomRanks);
    // RDKit❗✔️:         break;
    // RDKit❗✔️:       case RDKit::Atom::ChiralType::CHI_OCTAHEDRAL:
    // RDKit❗✔️:         embedOctahedral(mol, atm, efrags, atomRanks);
    // RDKit❗✔️:         break;
    // RDKit❗✔️:       default:
    // RDKit❗✔️:         break;
    // RDKit❗✔️:     }
    // RDKit❗✔️:   }
    // RDKit❗✔️: }
    // Behavior: source-order dispatch; incomplete geometry remains partial.
    // Complexity: one atom-rank pass and one bounded sort per chiral centre.
    let ranks: Vec<_> = (0..topology.atoms.len())
        .map(|index| atom_depict_rank(topology, index))
        .collect::<Result<_, _>>()?;
    let mut fragments = Vec::new();
    for (centre, atom) in topology.atoms.iter().enumerate() {
        let points = match atom.chiral_tag() {
            ChiralTag::SquarePlanar => square_planar_points(topology, centre, &ranks)?,
            ChiralTag::TrigonalBipyramidal => {
                trigonal_bipyramidal_points(topology, centre, &ranks)?
            }
            ChiralTag::Octahedral => octahedral_points(topology, centre, &ranks),
            _ => continue,
        };
        fragments.push(EmbeddedFrag::from_coord_map(topology, rings, &points)?);
    }
    Ok(fragments)
}
