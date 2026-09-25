//! Detached ligand relationships for non-tetrahedral stereocentres.

use cosmolkit_model::{ChiralTag, TopologyBlock};

const SQUAREPLANAR_ACROSS: [[u8; 4]; 4] = [[4, 4, 4, 4], [2, 3, 0, 1], [1, 0, 3, 2], [3, 2, 1, 0]];

const TRIGONALBIPYRAMIDAL_ACROSS: [[u8; 5]; 21] = [
    [5, 5, 5, 5, 5],
    [4, 5, 5, 5, 0],
    [4, 5, 5, 5, 0],
    [3, 5, 5, 0, 5],
    [3, 5, 5, 0, 5],
    [2, 5, 0, 5, 5],
    [2, 5, 0, 5, 5],
    [1, 0, 5, 5, 5],
    [1, 0, 5, 5, 5],
    [5, 4, 5, 5, 1],
    [5, 3, 5, 1, 5],
    [5, 4, 5, 5, 1],
    [5, 3, 5, 1, 5],
    [5, 2, 1, 5, 5],
    [5, 2, 1, 5, 5],
    [5, 5, 4, 5, 2],
    [5, 5, 3, 2, 5],
    [5, 5, 5, 4, 3],
    [5, 5, 5, 4, 3],
    [5, 5, 3, 2, 5],
    [5, 5, 4, 5, 2],
];

const OCTAHEDRAL_ACROSS: [[u8; 6]; 31] = [
    [6, 6, 6, 6, 6, 6],
    [5, 3, 4, 1, 2, 0],
    [5, 3, 4, 1, 2, 0],
    [4, 3, 5, 1, 0, 2],
    [5, 4, 3, 2, 1, 0],
    [4, 5, 3, 2, 0, 1],
    [3, 4, 5, 0, 1, 2],
    [3, 5, 4, 0, 2, 1],
    [5, 2, 1, 4, 3, 0],
    [4, 2, 1, 5, 0, 3],
    [5, 2, 1, 4, 3, 0],
    [4, 2, 1, 5, 0, 3],
    [3, 2, 1, 0, 5, 4],
    [3, 2, 1, 0, 5, 4],
    [5, 4, 3, 2, 1, 0],
    [4, 5, 3, 2, 0, 1],
    [4, 3, 5, 1, 0, 2],
    [3, 5, 4, 0, 2, 1],
    [3, 4, 5, 0, 1, 2],
    [2, 4, 0, 5, 1, 3],
    [2, 5, 0, 4, 3, 1],
    [2, 3, 0, 1, 5, 4],
    [2, 3, 0, 1, 5, 4],
    [2, 5, 0, 4, 3, 1],
    [2, 4, 0, 5, 1, 3],
    [1, 0, 4, 5, 2, 3],
    [1, 0, 5, 4, 3, 2],
    [1, 0, 3, 2, 5, 4],
    [1, 0, 3, 2, 5, 4],
    [1, 0, 5, 4, 3, 2],
    [1, 0, 4, 5, 2, 3],
];

const TRIGONALBIPYRAMIDAL_AXIAL: [[u8; 2]; 21] = [
    [5, 5],
    [0, 4],
    [0, 4],
    [0, 3],
    [0, 3],
    [0, 2],
    [0, 2],
    [0, 1],
    [0, 1],
    [1, 4],
    [1, 4],
    [1, 3],
    [1, 3],
    [1, 2],
    [1, 2],
    [2, 4],
    [2, 3],
    [3, 4],
    [3, 4],
    [2, 3],
    [2, 4],
];

/// The source's bond-order-relative across ligand, if present.
pub fn non_tetrahedral_across_ligand(
    topology: &TopologyBlock,
    centre: usize,
    ligand: usize,
) -> Option<usize> {
    // RDKit❗✔️: Bond *getChiralAcrossBond(const Atom *cen, const Bond *qry) {
    // RDKit❗✔️:   Atom::ChiralType tag = cen->getChiralTag();
    // RDKit❗✔️:   unsigned int perm = 0;
    // RDKit❗✔️:   cen->getPropIfPresent(common_properties::_chiralPermutation, perm);
    // RDKit❗✔️:   if (!perm) {
    // RDKit❗✔️:     return nullptr;
    // RDKit❗✔️:   }
    // RDKit❗✔️:   auto &mol = cen->getOwningMol();
    // RDKit❗✔️:   unsigned int count = 0;
    // RDKit❗✔️:   Bond *ref[6];
    // RDKit❗✔️:   int found = -1;
    // RDKit❗✔️:   unsigned int ref_max = getMaxNbors(tag);
    // RDKit❗✔️:   for (auto bnd : mol.atomBonds(cen)) {
    // RDKit❗✔️:     if (count == ref_max) {
    // RDKit❗✔️:       return nullptr;
    // RDKit❗✔️:     }
    // RDKit❗✔️:     ref[count] = bnd;
    // RDKit❗✔️:     if (bnd == qry) {
    // RDKit❗✔️:       found = count;
    // RDKit❗✔️:     }
    // RDKit❗✔️:     count++;
    // RDKit❗✔️:   }
    // RDKit❗✔️:   if (found >= 0) {
    // RDKit❗✔️:     switch (tag) {
    // RDKit❗✔️:       case Atom::ChiralType::CHI_SQUAREPLANAR:
    // RDKit❗✔️:         if (perm <= 3) {
    // RDKit❗✔️:           found = squareplanar_across[perm][found];
    // RDKit❗✔️:         } else {
    // RDKit❗✔️:           found = 4;
    // RDKit❗✔️:         }
    // RDKit❗✔️:         break;
    // RDKit❗✔️:       case Atom::ChiralType::CHI_TRIGONALBIPYRAMIDAL:
    // RDKit❗✔️:         if (perm <= 20) {
    // RDKit❗✔️:           found = trigonalbipyramidal_across[perm][found];
    // RDKit❗✔️:         } else {
    // RDKit❗✔️:           found = 5;
    // RDKit❗✔️:         }
    // RDKit❗✔️:         break;
    // RDKit❗✔️:       case Atom::ChiralType::CHI_OCTAHEDRAL:
    // RDKit❗✔️:         if (perm <= 30) {
    // RDKit❗✔️:           found = octahedral_across[perm][found];
    // RDKit❗✔️:         } else {
    // RDKit❗✔️:           found = 6;
    // RDKit❗✔️:         }
    // RDKit❗✔️:         break;
    // RDKit❗✔️:       default:
    // RDKit❗✔️:         return nullptr;
    // RDKit❗✔️:     }
    // RDKit❗✔️:     if (static_cast<unsigned int>(found) < count) {
    // RDKit❗✔️:       return ref[found];
    // RDKit❗✔️:     }
    // RDKit❗✔️:   }
    // RDKit❗✔️:   return nullptr;
    // RDKit❗✔️: }
    // Behavior: the tables above are copied from NontetrahedralStereo.cpp;
    // incomplete ligand sets may not contain the indicated across partner.
    // Complexity: one adjacency scan, matching the source's atomBonds walk.
    let atom = topology.atoms.get(centre)?;
    let permutation = usize::try_from(atom.chiral_permutation()?).ok()?;
    if permutation == 0 {
        return None;
    }
    let neighbors = topology.adjacency.neighbors_of(centre);
    let max_neighbors = match atom.chiral_tag() {
        ChiralTag::SquarePlanar => 4,
        ChiralTag::TrigonalBipyramidal => 5,
        ChiralTag::Octahedral => 6,
        _ => return None,
    };
    if neighbors.len() > max_neighbors {
        return None;
    }
    let found = neighbors
        .iter()
        .position(|neighbor| neighbor.atom_index == ligand)?;
    let across = match atom.chiral_tag() {
        ChiralTag::SquarePlanar => SQUAREPLANAR_ACROSS.get(permutation)?.get(found),
        ChiralTag::TrigonalBipyramidal => TRIGONALBIPYRAMIDAL_ACROSS.get(permutation)?.get(found),
        ChiralTag::Octahedral => OCTAHEDRAL_ACROSS.get(permutation)?.get(found),
        _ => None,
    }?;
    neighbors
        .get(usize::from(*across))
        .map(|neighbor| neighbor.atom_index)
}

/// The source's positive or negative axial ligand index, if present.
pub fn trigonal_bipyramidal_axial_ligand(
    topology: &TopologyBlock,
    centre: usize,
    axial: i32,
) -> Option<usize> {
    // RDKit❗✔️: Bond *getTrigonalBipyramidalAxialBond(const Atom *cen, int axial) {
    // RDKit❗✔️:   if (cen->getChiralTag() != RDKit::Atom::ChiralType::CHI_TRIGONALBIPYRAMIDAL ||
    // RDKit❗✔️:       cen->getDegree() > 5) {
    // RDKit❗✔️:     return nullptr;
    // RDKit❗✔️:   }
    // RDKit❗✔️:   unsigned int perm = 0;
    // RDKit❗✔️:   cen->getPropIfPresent(RDKit::common_properties::_chiralPermutation, perm);
    // RDKit❗✔️:   if (perm == 0 || perm > 20) {
    // RDKit❗✔️:     return nullptr;
    // RDKit❗✔️:   }
    // RDKit❗✔️:   unsigned int idx = (axial != -1) ? trigonalbipyramidal_axial[perm][0]
    // RDKit❗✔️:                                    : trigonalbipyramidal_axial[perm][1];
    // RDKit❗✔️:   unsigned int count = 0;
    // RDKit❗✔️:   for (const auto bnd : cen->getOwningMol().atomBonds(cen)) {
    // RDKit❗✔️:     if (count == idx) {
    // RDKit❗✔️:       return bnd;
    // RDKit❗✔️:     }
    // RDKit❗✔️:     count++;
    // RDKit❗✔️:   }
    // RDKit❗✔️:   return nullptr;
    // RDKit❗✔️: }
    // Behavior: keeps source permutation table and bond insertion ordering.
    // Complexity: indexed neighbor lookup replaces a linear source scan.
    let atom = topology.atoms.get(centre)?;
    let neighbors = topology.adjacency.neighbors_of(centre);
    if atom.chiral_tag() != ChiralTag::TrigonalBipyramidal || neighbors.len() > 5 {
        return None;
    }
    let permutation = usize::try_from(atom.chiral_permutation()?).ok()?;
    let pair = TRIGONALBIPYRAMIDAL_AXIAL.get(permutation)?;
    if permutation == 0 {
        return None;
    }
    let index = if axial == -1 { pair[1] } else { pair[0] };
    neighbors
        .get(usize::from(index))
        .map(|neighbor| neighbor.atom_index)
}

/// Source ideal ligand angle in degrees; zero denotes a non-applicable tag.
pub fn non_tetrahedral_ideal_angle(
    topology: &TopologyBlock,
    centre: usize,
    first: usize,
    second: usize,
) -> f64 {
    // RDKit❗✔️: double getIdealAngleBetweenLigands(const Atom *cen, const Atom *lig1,
    // RDKit❗✔️:                                    const Atom *lig2) {
    // RDKit❗✔️:   auto tag = cen->getChiralTag();
    // RDKit❗✔️:   switch (tag) {
    // RDKit❗✔️:     case Atom::ChiralType::CHI_SQUAREPLANAR:
    // RDKit❗✔️:     case Atom::ChiralType::CHI_OCTAHEDRAL:
    // RDKit❗✔️:       return getChiralAcrossAtom(cen, lig1) == lig2 ? 180 : 90;
    // RDKit❗✔️:       break;
    // RDKit❗✔️:     case Atom::ChiralType::CHI_TRIGONALBIPYRAMIDAL:
    // RDKit❗✔️:       if (getChiralAcrossAtom(cen, lig1) == lig2) {
    // RDKit❗✔️:         return 180;
    // RDKit❗✔️:       } else if (isTrigonalBipyramidalAxialAtom(cen, lig1) ||
    // RDKit❗✔️:                  isTrigonalBipyramidalAxialAtom(cen, lig2)) {
    // RDKit❗✔️:         return 90;
    // RDKit❗✔️:       } else {
    // RDKit❗✔️:         return 120;
    // RDKit❗✔️:       }
    // RDKit❗✔️:       break;
    // RDKit❗✔️:     default:
    // RDKit❗✔️:       return 0;
    // RDKit❗✔️:   }
    // RDKit❗✔️: }
    // Behavior: even a missing across ligand compares as false to a real ligand.
    // Complexity: bounded (six-neighbor) adjacency scans, as in source.
    let tag = topology.atoms.get(centre).map(|atom| atom.chiral_tag());
    let across = non_tetrahedral_across_ligand(topology, centre, first);
    match tag {
        Some(ChiralTag::SquarePlanar | ChiralTag::Octahedral) => {
            if across == Some(second) {
                180.0
            } else {
                90.0
            }
        }
        Some(ChiralTag::TrigonalBipyramidal) => {
            if across == Some(second) {
                180.0
            } else if trigonal_bipyramidal_axial_ligand(topology, centre, 1) == Some(first)
                || trigonal_bipyramidal_axial_ligand(topology, centre, -1) == Some(first)
                || trigonal_bipyramidal_axial_ligand(topology, centre, 1) == Some(second)
                || trigonal_bipyramidal_axial_ligand(topology, centre, -1) == Some(second)
            {
                90.0
            } else {
                120.0
            }
        }
        _ => 0.0,
    }
}
