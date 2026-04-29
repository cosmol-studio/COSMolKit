use std::collections::{HashMap, VecDeque};

use crate::{Bond, BondOrder, ChiralTag, Molecule, ValenceModel, assign_valence};

/// A ligand participating in a tetrahedral stereocenter.
#[derive(Debug, Copy, Clone, PartialEq, Eq)]
pub enum LigandRef {
    /// A graph atom by molecule atom index.
    Atom(usize),
    /// A hydrogen carried as an atom hydrogen count rather than a graph atom.
    ImplicitH,
}

/// Ordered tetrahedral stereochemistry for one atom center.
///
/// Specification: `tetrahedral_stereo_representation.md` at repository root.
///
/// For RDKit-style tetrahedral tags, the ordering is derived from the atom's
/// incoming bond order and chiral tag.
#[derive(Debug, Copy, Clone, PartialEq, Eq)]
pub struct TetrahedralStereo {
    pub center: usize,
    pub ligands: [LigandRef; 4],
}

/// Hidden helper struct exposing RDKit legacy stereochemistry perception
/// fields for parity tests.
#[doc(hidden)]
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct LegacyStereoAtomProps {
    pub cip_code: Option<String>,
    pub cip_rank: Option<i64>,
    pub chirality_possible: bool,
}

fn bond_affects_atom_chirality(bond: &Bond, atom_index: usize) -> bool {
    !matches!(bond.order, BondOrder::Null)
        && !(matches!(bond.order, BondOrder::Dative) && bond.begin_atom == atom_index)
}

fn atom_nonzero_degree(mol: &Molecule, atom_index: usize) -> usize {
    mol.bonds
        .iter()
        .filter(|bond| {
            (bond.begin_atom == atom_index || bond.end_atom == atom_index)
                && bond_affects_atom_chirality(bond, atom_index)
        })
        .count()
}

fn has_protium_neighbor(mol: &Molecule, atom_index: usize) -> bool {
    mol.bonds.iter().any(|bond| {
        let nbr = if bond.begin_atom == atom_index {
            Some(bond.end_atom)
        } else if bond.end_atom == atom_index {
            Some(bond.begin_atom)
        } else {
            None
        };
        nbr.is_some_and(|nbr_idx| {
            let atom = &mol.atoms[nbr_idx];
            atom.atomic_num == 1 && atom.isotope.is_none()
        })
    })
}

fn bond_is_conjugated_for_stereo(mol: &Molecule, bond_index: usize) -> bool {
    let bond = &mol.bonds[bond_index];
    if !matches!(
        bond.order,
        BondOrder::Single | BondOrder::Double | BondOrder::Aromatic
    ) {
        return false;
    }
    let endpoints = [bond.begin_atom, bond.end_atom];
    endpoints.into_iter().any(|center| {
        let adjacent: Vec<usize> = mol
            .bonds
            .iter()
            .enumerate()
            .filter_map(|(idx, other)| {
                if idx == bond_index {
                    return None;
                }
                if other.begin_atom == center || other.end_atom == center {
                    Some(idx)
                } else {
                    None
                }
            })
            .collect();
        adjacent.iter().any(|&a| {
            adjacent.iter().any(|&b| {
                a != b
                    && matches!(
                        (mol.bonds[a].order, mol.bonds[b].order),
                        (BondOrder::Single, BondOrder::Double)
                            | (BondOrder::Double, BondOrder::Single)
                            | (BondOrder::Aromatic, BondOrder::Single)
                            | (BondOrder::Single, BondOrder::Aromatic)
                            | (BondOrder::Aromatic, BondOrder::Double)
                            | (BondOrder::Double, BondOrder::Aromatic)
                            | (BondOrder::Aromatic, BondOrder::Aromatic)
                    )
            })
        })
    })
}

fn atom_has_conjugated_bond(mol: &Molecule, atom_index: usize) -> bool {
    mol.bonds.iter().enumerate().any(|(bond_index, bond)| {
        (bond.begin_atom == atom_index || bond.end_atom == atom_index)
            && bond_is_conjugated_for_stereo(mol, bond_index)
    })
}

fn count_swaps_to_interconvert(probe: &[usize], reference: &[usize]) -> usize {
    let mut probe = probe.to_vec();
    let mut positions: HashMap<usize, usize> =
        probe.iter().enumerate().map(|(i, &v)| (v, i)).collect();
    let mut swaps = 0usize;
    for (i, &target) in reference.iter().enumerate() {
        if probe[i] == target {
            continue;
        }
        let j = *positions
            .get(&target)
            .expect("reference element missing from probe permutation");
        let current = probe[i];
        probe.swap(i, j);
        positions.insert(current, j);
        positions.insert(target, i);
        swaps += 1;
    }
    swaps
}

fn perturbation_order_for_atom(mol: &Molecule, atom_index: usize, probe: &[usize]) -> usize {
    let reference: Vec<usize> = mol
        .bonds
        .iter()
        .filter(|bond| bond.begin_atom == atom_index || bond.end_atom == atom_index)
        .map(|bond| bond.index)
        .collect();
    count_swaps_to_interconvert(probe, &reference)
}

fn min_cycle_size_for_bond(mol: &Molecule, bond_index: usize) -> Option<usize> {
    let bond = &mol.bonds[bond_index];
    let start = bond.begin_atom;
    let goal = bond.end_atom;
    let mut queue = VecDeque::from([(start, 0usize)]);
    let mut seen = vec![false; mol.atoms.len()];
    seen[start] = true;
    while let Some((atom, dist)) = queue.pop_front() {
        for (idx, edge) in mol.bonds.iter().enumerate() {
            if idx == bond_index {
                continue;
            }
            let next = if edge.begin_atom == atom {
                Some(edge.end_atom)
            } else if edge.end_atom == atom {
                Some(edge.begin_atom)
            } else {
                None
            };
            let Some(next_atom) = next else {
                continue;
            };
            if next_atom == goal {
                return Some(dist + 2);
            }
            if !seen[next_atom] {
                seen[next_atom] = true;
                queue.push_back((next_atom, dist + 1));
            }
        }
    }
    None
}

fn should_detect_double_bond_stereo(mol: &Molecule, bond_index: usize) -> bool {
    min_cycle_size_for_bond(mol, bond_index).is_none_or(|size| size >= 8)
}

fn is_atom_potential_chiral_center(
    mol: &Molecule,
    atom_index: usize,
    ranks: &[i64],
    explicit_valence: &[u8],
    implicit_hydrogens: &[u8],
) -> (bool, bool, Vec<(i64, usize)>) {
    let atom = &mol.atoms[atom_index];
    let nz_degree = atom_nonzero_degree(mol, atom_index);
    let graph_h_neighbors = mol.bonds.iter().filter(|bond| {
        let nbr = if bond.begin_atom == atom_index {
            Some(bond.end_atom)
        } else if bond.end_atom == atom_index {
            Some(bond.begin_atom)
        } else {
            None
        };
        nbr.is_some_and(|nbr_idx| mol.atoms[nbr_idx].atomic_num == 1)
    });
    let total_num_hs = atom.explicit_hydrogens as usize
        + implicit_hydrogens[atom_index] as usize
        + graph_h_neighbors.count();
    let tnz_degree = nz_degree + total_num_hs;
    let mut legal_center = true;
    let mut has_dupes = false;
    let mut nbrs = Vec::new();

    if tnz_degree > 4 {
        legal_center = false;
    } else if tnz_degree < 3 {
        legal_center = false;
    } else if nz_degree < 3 && !matches!(atom.atomic_num, 15 | 33) {
        legal_center = false;
    } else if nz_degree == 3 {
        if total_num_hs == 1 {
            if has_protium_neighbor(mol, atom_index) {
                legal_center = false;
            }
        } else {
            legal_center = false;
            match atom.atomic_num {
                7 => {
                    let in_three_ring = mol.bonds.iter().enumerate().any(|(bond_index, bond)| {
                        (bond.begin_atom == atom_index || bond.end_atom == atom_index)
                            && min_cycle_size_for_bond(mol, bond_index) == Some(3)
                    });
                    if !atom_has_conjugated_bond(mol, atom_index) && in_three_ring {
                        legal_center = true;
                    }
                }
                15 | 33 => {
                    legal_center = true;
                }
                16 | 34 => {
                    let ev = explicit_valence[atom_index] as i32;
                    if ev == 4 || (ev == 3 && atom.formal_charge == 1) {
                        legal_center = true;
                    }
                }
                _ => {}
            }
        }
    }

    if legal_center && !ranks.is_empty() {
        let mut seen_ranks = HashMap::<i64, ()>::new();
        for bond in &mol.bonds {
            let other_idx = if bond.begin_atom == atom_index {
                Some(bond.end_atom)
            } else if bond.end_atom == atom_index {
                Some(bond.begin_atom)
            } else {
                None
            };
            let Some(other_idx) = other_idx else {
                continue;
            };
            nbrs.push((ranks[other_idx], bond.index));
            if !bond_affects_atom_chirality(bond, atom_index) {
                continue;
            }
            if seen_ranks.insert(ranks[other_idx], ()).is_some() {
                has_dupes = true;
                break;
            }
        }
    }

    (legal_center, has_dupes, nbrs)
}

impl Molecule {
    /// Return ordered tetrahedral stereocenters derived from RDKit-style atom
    /// chiral tags and molecule bond order.
    #[must_use]
    pub fn tetrahedral_stereo(&self) -> Vec<TetrahedralStereo> {
        let mut out = Vec::new();
        for atom in &self.atoms {
            let tag = atom.chiral_tag;
            if !matches!(tag, ChiralTag::TetrahedralCcw | ChiralTag::TetrahedralCw) {
                continue;
            }

            let mut ligands = Vec::with_capacity(4);
            for bond in &self.bonds {
                if bond.begin_atom == atom.index {
                    ligands.push(LigandRef::Atom(bond.end_atom));
                } else if bond.end_atom == atom.index {
                    ligands.push(LigandRef::Atom(bond.begin_atom));
                }
            }

            if ligands.len() == 3 {
                ligands.push(LigandRef::ImplicitH);
            }
            let Ok(mut ligands) = <[LigandRef; 4]>::try_from(ligands) else {
                continue;
            };

            // RDKit DistGeom treats CHI_TETRAHEDRAL_CCW as positive volume
            // for atom bond order, and CHI_TETRAHEDRAL_CW as negative. COSMolKit
            // stores ordered ligands with positive orientation, so CW is an odd
            // permutation of the RDKit atom-bond order.
            if matches!(tag, ChiralTag::TetrahedralCw) {
                ligands.swap(0, 1);
            }

            out.push(TetrahedralStereo {
                center: atom.index,
                ligands,
            });
        }
        out
    }

    /// Hidden helper that mirrors the atom-property side of RDKit legacy
    /// stereochemistry perception closely enough for parity testing.
    #[doc(hidden)]
    #[must_use]
    pub fn rdkit_legacy_stereo_atom_props(
        &self,
        flag_possible_stereo_centers: bool,
    ) -> Vec<LegacyStereoAtomProps> {
        let Ok(assignment) = assign_valence(self, ValenceModel::RdkitLike) else {
            return vec![
                LegacyStereoAtomProps {
                    cip_code: None,
                    cip_rank: None,
                    chirality_possible: false,
                };
                self.atoms.len()
            ];
        };

        let mut has_stereo_atoms = false;
        let mut has_potential_stereo_atoms = false;
        for atom in &self.atoms {
            if !has_stereo_atoms && !matches!(atom.chiral_tag, ChiralTag::Unspecified) {
                has_stereo_atoms = true;
            } else if !has_potential_stereo_atoms {
                has_potential_stereo_atoms = is_atom_potential_chiral_center(
                    self,
                    atom.index,
                    &[],
                    &assignment.explicit_valence,
                    &assignment.implicit_hydrogens,
                )
                .0;
            }
        }

        let mut has_stereo_bonds = false;
        let mut has_potential_stereo_bonds = false;
        for (bond_index, bond) in self.bonds.iter().enumerate() {
            if !matches!(bond.order, BondOrder::Double) {
                continue;
            }
            let begin = bond.begin_atom;
            let end = bond.end_atom;
            let is_specified = self.bonds.iter().any(|other| {
                (other.begin_atom == begin
                    || other.end_atom == begin
                    || other.begin_atom == end
                    || other.end_atom == end)
                    && matches!(
                        other.direction,
                        crate::BondDirection::EndUpRight | crate::BondDirection::EndDownRight
                    )
            });
            if is_specified {
                has_stereo_bonds = true;
            } else if !has_potential_stereo_bonds
                && should_detect_double_bond_stereo(self, bond_index)
            {
                has_potential_stereo_bonds = true;
            }
            if has_stereo_bonds && has_potential_stereo_bonds {
                break;
            }
        }

        let mut keep_going = has_stereo_atoms || has_stereo_bonds;
        if !keep_going {
            keep_going = flag_possible_stereo_centers
                && (has_potential_stereo_atoms || has_potential_stereo_bonds);
        }
        if !keep_going {
            return vec![
                LegacyStereoAtomProps {
                    cip_code: None,
                    cip_rank: None,
                    chirality_possible: false,
                };
                self.atoms.len()
            ];
        }

        let ranks = crate::io::molblock::rdkit_cip_ranks_for_depict(self);
        let mut out: Vec<LegacyStereoAtomProps> = ranks
            .iter()
            .map(|&rank| LegacyStereoAtomProps {
                cip_code: None,
                cip_rank: Some(rank),
                chirality_possible: false,
            })
            .collect();

        for atom in &self.atoms {
            let tag = atom.chiral_tag;
            if !flag_possible_stereo_centers && matches!(tag, ChiralTag::Unspecified) {
                continue;
            }

            let (legal_center, has_dupes, mut nbrs) = is_atom_potential_chiral_center(
                self,
                atom.index,
                &ranks,
                &assignment.explicit_valence,
                &assignment.implicit_hydrogens,
            );
            if legal_center && !has_dupes && flag_possible_stereo_centers {
                out[atom.index].chirality_possible = true;
            }
            if !legal_center || has_dupes || matches!(tag, ChiralTag::Unspecified) {
                continue;
            }

            nbrs.sort();
            let probe: Vec<usize> = nbrs.into_iter().map(|(_, bond_index)| bond_index).collect();
            let mut swaps = perturbation_order_for_atom(self, atom.index, &probe);
            let total_num_hs = atom.explicit_hydrogens as usize
                + assignment.implicit_hydrogens[atom.index] as usize
                + self
                    .bonds
                    .iter()
                    .filter(|bond| {
                        let nbr = if bond.begin_atom == atom.index {
                            Some(bond.end_atom)
                        } else if bond.end_atom == atom.index {
                            Some(bond.begin_atom)
                        } else {
                            None
                        };
                        nbr.is_some_and(|nbr_idx| self.atoms[nbr_idx].atomic_num == 1)
                    })
                    .count();
            if probe.len() == 3 && total_num_hs == 1 {
                swaps += 1;
            }
            let mut final_tag = tag;
            if swaps % 2 == 1 {
                final_tag = match final_tag {
                    ChiralTag::TetrahedralCcw => ChiralTag::TetrahedralCw,
                    ChiralTag::TetrahedralCw => ChiralTag::TetrahedralCcw,
                    ChiralTag::Unspecified => ChiralTag::Unspecified,
                };
            }
            out[atom.index].cip_code = Some(
                if matches!(final_tag, ChiralTag::TetrahedralCcw) {
                    "S"
                } else {
                    "R"
                }
                .to_owned(),
            );
        }
        out
    }

    /// Hidden helper mirroring RDKit's legacy small-ring gate for double-bond
    /// stereochemistry.
    #[doc(hidden)]
    #[must_use]
    pub fn rdkit_should_detect_double_bond_stereo(&self, bond_index: usize) -> bool {
        should_detect_double_bond_stereo(self, bond_index)
    }
}

#[cfg(test)]
mod tests {
    use super::{LigandRef, TetrahedralStereo};
    use crate::Molecule;

    #[test]
    fn tetrahedral_stereo_is_derived_from_chiral_tags() {
        let ccw = Molecule::from_smiles("F[C@](Cl)(Br)I").expect("parse chiral smiles");
        let cw = Molecule::from_smiles("F[C@@](Cl)(Br)I").expect("parse chiral smiles");

        let ccw_stereo = ccw.tetrahedral_stereo();
        let cw_stereo = cw.tetrahedral_stereo();

        assert_eq!(ccw_stereo.len(), 1);
        assert_eq!(cw_stereo.len(), 1);
        assert_eq!(ccw_stereo[0].center, 1);
        assert_eq!(cw_stereo[0].center, 1);
        assert_eq!(
            ccw_stereo[0].ligands,
            [
                LigandRef::Atom(0),
                LigandRef::Atom(2),
                LigandRef::Atom(3),
                LigandRef::Atom(4)
            ]
        );
        assert_eq!(
            cw_stereo[0].ligands,
            [
                LigandRef::Atom(2),
                LigandRef::Atom(0),
                LigandRef::Atom(3),
                LigandRef::Atom(4)
            ]
        );
    }

    #[test]
    fn tetrahedral_stereo_places_implicit_hydrogen_as_fourth_ligand() {
        let mol = Molecule::from_smiles("[13CH3:7][C@H](F)Cl").expect("parse chiral smiles");
        let stereo = mol.tetrahedral_stereo();

        assert_eq!(
            stereo,
            vec![TetrahedralStereo {
                center: 1,
                ligands: [
                    LigandRef::Atom(0),
                    LigandRef::Atom(2),
                    LigandRef::Atom(3),
                    LigandRef::ImplicitH
                ],
            }]
        );
    }
}
